function PAT_breast

load PAT/PAT_breast_boundary.mat
load PAT/PAT_setup.mat
hmin = .25;
% hmin = .2;
% hmin = .175;

tb = linspace(min(t),max(t)-hmin,ceil((max(t)-min(t))/hmin));
xB = interp1(t,xb,tb);
yB = interp1(t,yb,tb);
tb = [tb tb(1)]; 
xB = [xB xB(1)]; 
yB = [yB yB(1)]; 

node = [xB(:) yB(:)];
hdata.fun = @(x,y) hmin*ones(size(x));
hdata.hmax = .025*hmin;
[p,tri] = mesh2d(node,[],hdata);

d1 = size(speed,1); d2 = size(speed,2);
scale = .2*.001; % .2mm per pixel * millimeters to meters factor
a = d1*scale;
b = d2*scale;
[xm ym] = meshgrid(linspace(0,a,d1),linspace(0,b,d2));
pcolor(xm,ym,speed');shading interp
hold on

plot(xB,yB,'ko--','markersize',16);
plot(p(:,1),p(:,2),'ro','markersize',8)
hold on;
triplot(tri,p(:,1),p(:,2),'k','linewidth',2)
axis on

%% produce model of c2

Globals2D
N = 2;
FinalTime = 1e-4;

% make mesh
VX = p(:,1); VX = VX(:)';
VY = p(:,2); VY = VY(:)'; 
EToV = tri;
K = size(EToV,1);
StartUp2D

% move extra vertices to boundary
fids = {[1 2],[2 3],[3 1]};
for e = 1:K
    ev = EToV(e,:);
    for f = 1:Nfaces
        if EToE(e,f)==e
            vxb = VX(ev(fids{f}));
            vyb = VY(ev(fids{f}));
            for i = 1:2
                [val id] = min(abs(xB - vxb(i)) + abs(yB - vyb(i)));
                if val > 1e-8                                       
                    newt = atan2(vyb(i)-mean(yb),vxb(i)-mean(xb));
                    newx = interp1(t,xb,newt);
                    newy = interp1(t,yb,newt);
                    VX(ev(fids{f}(i))) = newx;
                    VY(ev(fids{f}(i))) = newy;                    
                end                    
            end
        end
    end
end
StartUp2D

% set up cubature/plotting
global Pq Vq cq rho
[rq sq wq] = Cubature2D(2*N);
Vq = Vandermonde2D(N,rq,sq)/V;
Pq = V*V' * Vq'*diag(wq);
xq = Vq*x;
yq = Vq*y;
Jq = Vq*J;

[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x;
yp = Vp*y;

% smooth velocity model
for iter = 1:5
    savg = speed;
    for i = 2:d1-1
        for j = 2:d2-1
            sij = speed(i,j) + ...
                speed(i+1,j) + speed(i-1,j) + speed(i,j+1) + speed(i,j-1) + ...
                speed(i+1,j+1) + speed(i-1,j+1) + speed(i+1,j-1) + speed(i-1,j-1);
            savg(i,j) = sij/9;
        end
    end
    speed = savg;
end

cq = interp2(xm,ym,speed',xq,yq);
rho = interp2(xm,ym,density',xq,yq);

% clf
% vv = cq;
% color_line3(xq,yq,vv,vv,'.')
% % vv = Vp*c2; color_line3(xp,yp,vv,vv,'.')
% hold on
% axis on
% c2 = Pq*cq;
% plot3(x(Fmask(:),:),y(Fmask(:),:),c2(Fmask(:),:)*1.05,'k-','linewidth',2)
% colorbar
% return


%% init cond

% smooth absorption initial cond
for iter = 1:5
    savg = absorption;
    for i = 2:d1-1
        for j = 2:d2-1
            sij = absorption(i,j) + ...
                absorption(i+1,j) + absorption(i-1,j) + absorption(i,j+1) + absorption(i,j-1) + ...
                absorption(i+1,j+1) + absorption(i-1,j+1) + absorption(i+1,j-1) + absorption(i-1,j-1);
            savg(i,j) = sij/9;
        end
    end
    absorption = savg;
end

pex = @(x,y) interp2(xm,ym,absorption',x,y);

% projection
p = Pq*pex(xq,yq);
u = zeros(Np, K);
v = zeros(Np, K);

% keyboard
% vv = Vp*p; color_line3(xp,yp,vv,vv,'.'); axis equal; axis tight; colorbar
% return

%% check eigs
if 0 & 3*Np*K < 4000
    U = zeros(Np*K,3); A = zeros(3*Np*K);
    for i = 1:3*Np*K
        U(i) = 1;
        p = reshape(U(:,1),Np,K); u = reshape(U(:,2),Np,K); v = reshape(U(:,3),Np,K);
        [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,0);
        A(:,i) = [rhsp(:); rhsu(:); rhsv(:)];
        U(i) = 0;
        if mod(i,round(3*Np*K/10))==0
            disp(sprintf('i = %d out of %d\n',i,3*Np*K))
        end
    end
    lam = eig(A);
    keyboard
end

%% setup timestepper

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
CNh = max(max(cq)*CN*max(Fscale(:)));
dt = 2/CNh;

Nstep = ceil(FinalTime/dt);
dt = FinalTime/Nstep;

tau = 1;

%% generate synthetic data

NstepRK = Nstep*5;
Ubc{1} = zeros(length(mapB(:)),NstepRK);
Ubc{2} = zeros(length(mapB(:)),NstepRK);

% outer time step loop
time = 0;
for i = 1:Nstep
    for INTRK = 1:5
        [rhsp, rhsu, rhsv,UbcI] = acousticsRHS2D(p,u,v,tau,'abc');
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        p = p+rk4b(INTRK)*resp;
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        
        % save traces
        id0 = ((i-1)*5+(INTRK-1)); % zero index
        Ubc{1}(:,NstepRK - id0) = UbcI{1};
        Ubc{2}(:,NstepRK - id0) = UbcI{2};
    end;
    
    if mod(i,10)==0
        vv = Vp*p;
        clf; color_line3(xp,yp,vv,vv,'.');
        axis equal; axis tight; colorbar
        title(sprintf('time = %f, max solution val = %f',time,max(abs(vv(:)))))
        drawnow
    end
    time = i * dt;
    if (mod(i,ceil(Nstep/10))==0)
        disp(sprintf('generating synthetic data: tstep %d out of %d\n',i,Nstep))
    end
end


disp('Synthetic data generated. ')


%% playback : reverse propagate to generate Rm from pressure boundary data

dt = -dt;
p = zeros(Np,K); u = p; v = p;
resp = zeros(Np,K); resu = resp; resv = resp;

for i = 1:Nstep
    for INTRK = 1:5
        id0 = ((i-1)*5+(INTRK-1)); % 1-index
        UbcI{1} = Ubc{1}(:,id0+1);
        UbcI{2} = Ubc{2}(:,id0+1);
        
        [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,sign(dt)*tau,'rd',UbcI);
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        p = p + rk4b(INTRK)*resp;
        u = u + rk4b(INTRK)*resu;
        v = v + rk4b(INTRK)*resv;
    end;
    
    if (mod(i,ceil(Nstep/5))==0)
        disp(sprintf('first backwards prop: tstep %d out of %d',i,Nstep))
    end
    if mod(i,10)==0
        vv = Vp*p;
        clf; color_line3(xp,yp,vv,vv,'.');
        axis equal; axis tight; colorbar
        title(sprintf('first reversal time = %f, max solution val = %f',time,max(abs(vv(:)))))
        drawnow
    end
end

% save reversed data
Rm = p;


%% iterate

% initalize pressure condition
p0 = Rm;

for iter = 1:5
    
    p = p0; u = zeros(Np,K); v = zeros(Np,K);
    resp = zeros(Np,K); resu = resp; resv = resp;
    
    % compute forward propagation
    dt = abs(dt);
    for i = 1:Nstep
        for INTRK = 1:5
            [rhsp, rhsu, rhsv, UbcI] = acousticsRHS2D(p,u,v,sign(dt)*tau,'abc');
            resp = rk4a(INTRK)*resp + dt*rhsp;
            resu = rk4a(INTRK)*resu + dt*rhsu;
            resv = rk4a(INTRK)*resv + dt*rhsv;
            p = p + rk4b(INTRK)*resp;
            u = u + rk4b(INTRK)*resu;
            v = v + rk4b(INTRK)*resv;
        end;
        
        if (mod(i,ceil(Nstep/5))==0)
            disp(sprintf('fwd prop: tstep %d out of %d',i,Nstep))
        end
        %         if mod(i,10)==0
        %             vv = Vp*p;
        %             clf; color_line3(xp,yp,vv,vv,'.');
        %             title(sprintf('forwards time = %f, max solution val = %f',time,max(abs(vv(:)))))
        %             axis equal; axis tight; colorbar
        %             drawnow
        %         end
    end
    
    %     % compute new auxiliary initial conditions
    %     p = 0*p;
    %     u = 0*u;
    %     v = 0*v;
    
    % compute backwards propagation
    dt = -dt;
    for i = 1:Nstep
        for INTRK = 1:5
            [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,sign(dt)*tau,'d');
            
            resp = rk4a(INTRK)*resp + dt*rhsp;
            resu = rk4a(INTRK)*resu + dt*rhsu;
            resv = rk4a(INTRK)*resv + dt*rhsv;
            p = p + rk4b(INTRK)*resp;
            u = u + rk4b(INTRK)*resu;
            v = v + rk4b(INTRK)*resv;
        end;
        
        if (mod(i,ceil(Nstep/5))==0)
            disp(sprintf('backwards prop: tstep %d out of %d',i,Nstep))
        end
        %         if mod(i,10)==0
        %             vv = Vp*p;
        %             clf; color_line3(xp,yp,vv,vv,'.');
        %             title(sprintf('backwards, time = %f, max solution val = %f',time,max(abs(vv(:)))))
        %             axis equal; axis tight; colorbar
        %             drawnow
        %         end
    end
    
    % update initial condition
    %     p0 = (p0-p) + Rm;
    p0 = p + Rm;
    
    % plot initial cond
    err = diag(wq)*(Vq*J).*(Vq*p0-pex(xq,yq)).^2;
    L2err(iter) = sqrt(sum(err(:)));
    uq = diag(wq)*(Vq*J).*pex(xq,yq).^2;
    Unorm = sqrt(sum(uq(:)));
    
    disp(['iter ' num2str(iter) ', L2 err in recon = ',num2str(L2err/Unorm)])
    %     pause
end

keyboard

figure
% vv = pex(xp,yp)-Vp*p0;
vv = Vp*p0;
% vv = pex(xp,yp);
clf; color_line3(xp,yp,vv,vv,'.');
axis equal; axis tight; colorbar
title(sprintf('iter = %d, reconstruction error = %f',iter,L2err(iter)))

return


function [rhsp, rhsu, rhsv, UbcO] = acousticsRHS2D(p,u,v,tau,bcOpt,Ubc)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

% evaluate upwind fluxes
ndotdU = nx.*du + ny.*dv;

if strcmp(bcOpt,'d') % Impose reflective boundary conditions (p+ = -p-)
    ndotdU(mapB) = 0;
    dp(mapB) = 2*(-p(vmapB));
elseif strcmp(bcOpt,'n') % free surface BCs
    ndotdU(mapB) = -2*(nx(mapB).*u(vmapM(mapB)) + ny(mapB).*v(vmapM(mapB)));
    dp(mapB) = 0;
elseif strcmp(bcOpt,'abc') % basic abcs
    un = nx(mapB).*u(vmapM(mapB)) + ny(mapB).*v(vmapM(mapB));
    dp(mapB) = -p(vmapB);
    ndotdU(mapB) = -un;
end

if strcmp(bcOpt,'rd')
    dp(mapB) = 2*(Ubc{1}-p(vmapB));
    ndotdU(mapB) = 0;
elseif strcmp(bcOpt,'rn')
    dp(mapB) = 0;
    ndotdU(mapB) = 2*(Ubc{2}-nx(mapB).*u(vmapM(mapB)) + ny(mapB).*v(vmapM(mapB)));
end
UbcO{1} = p(vmapB);
UbcO{2} = ndotdU(mapB);

fluxp =  tau*dp - ndotdU;
fluxu =  (tau*ndotdU - dp);

pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
rhsp =  -divU + LIFT*(Fscale.*fluxp)/2.0;
rhsu =  -dpdx + LIFT*(Fscale.*fluxu.*nx)/2.0;
rhsv =  -dpdy + LIFT*(Fscale.*fluxu.*ny)/2.0;

global Pq Vq cq rho
rhsp = Pq*(cq.*(Vq*rhsp));
rhsu = Pq*(rho.*(Vq*rhsu));
rhsv = Pq*(rho.*(Vq*rhsv));

return;

