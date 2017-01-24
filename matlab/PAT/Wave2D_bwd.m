function Wave2D

% clear all, clear
clear -global *

Globals2D

N = 5;
K1D = 8;
c_flag = 0;
cfun = @(x,y) ones(size(x));
% cfun = @(x,y) 1 + (x > 0);

% cfun = @(x,y) 1 + .25*sin(2*pi*x).*sin(2*pi*y); % smooth velocity
% cfun = @(x,y) (1 + .25*sin(2*pi*x).*sin(2*pi*y) + (y > 0)); % piecewise smooth velocity

% filename = 'Grid/Other/block2.neu';
% filename = 'Grid/Maxwell2D/Maxwell05.neu';
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
StartUp2D;

global Pq cq Vq

% PlotMesh2D; return

[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
Pq = V*V'*Vq'*diag(wq); % J's cancel out
Mref = inv(V*V');
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

%%
cq = cfun(xq,yq);

% clf; vv = Vp*Pq*cq; color_line3(xp,yp,vv,vv,'.');axis equal;axis tight;return

% FinalTime = 2/min(cq(:));
FinalTime = 2;

%% params setup

x0 = .1; y0 = 0;

% gaussian pulse
pex = @(x,y) exp(-5^2*((x-x0).^2 + (y-y0).^2));
% pex = @(x,y) cos(.5*pi*x).*cos(.5*pi*y);

% % smooth annulus
% sig = @(x) 1-1./(1+exp(-100*x));
% r2 = @(x,y) x.^2 + y.^2;
% pex = @(x,y) sig((r2(x,y)-.25).^2);

% % boxes overlapping
% a = 2/K1D;
% pex = @(x,y) (abs(x+a)<.5 & abs(y)<.5) + (abs(x)<.5 & abs(y-a)<.5) + (abs(x-a)<.5 & abs(y+a)<.5);

% projection
p = Pq*pex(xq,yq);
u = zeros(Np, K);
v = zeros(Np, K);

% vv = Vp*p; color_line3(xp,yp,vv,vv,'.'); axis equal; axis tight; return

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
CNh = max(CN*max(Fscale(:)));
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
    
%     if 1 && nargin==0 && mod(i,10)==0
%         vv = Vp*p;
%         clf; color_line3(xp,yp,vv,vv,'.');
%         axis equal; axis tight; colorbar        
%         title(sprintf('time = %f, max solution val = %f',time,max(abs(vv(:)))))
%         drawnow
%     end
    time = i * dt;
    if (mod(i,ceil(Nstep/10))==0)
        disp(sprintf('generating synthetic data: tstep %d out of %d\n',i,Nstep))
    end
end

disp('Synthetic data generated. ')
% pause

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
%     if mod(i,10)==0
%         vv = Vp*p;
%         clf; color_line3(xp,yp,vv,vv,'.');
%         axis equal; axis tight; colorbar        
%         title(sprintf('first reversal time = %f, max solution val = %f',time,max(abs(vv(:)))))
%         drawnow
%     end        
end

% save reversed data
Rm = p; 


%% iterate 

% initalize pressure condition
p0 = Rm;

for iter = 1:50
    
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
            
            % save traces
            id0 = ((i-1)*5+(INTRK-1)); % zero index
            Ubc{1}(:,NstepRK - id0) = UbcI{1};
            Ubc{2}(:,NstepRK - id0) = UbcI{2};

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
    
    % compute backwards propagation
    dt = -dt;
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
    p0 = (p0-p) + Rm; 
    
    % plot initial cond
    err = diag(wq)*(Vq*J).*(Vq*p0-pex(xq,yq)).^2;
    L2err(iter) = sqrt(sum(err(:)));   
    uq = diag(wq)*(Vq*J).*pex(xq,yq).^2;    
    Unorm = sqrt(sum(uq(:)));
%     figure
%     vv = pex(xp,yp)-Vp*p0;
%     clf; color_line3(xp,yp,vv,vv,'.');
%     axis equal; axis tight; colorbar
%     s = sprintf('iter = %d, reconstruction error = %f',iter,L2err(iter));
%     title(s)
    disp(['L2 err in recon = ', num2str(L2err/Unorm)])
%     pause
end

keyboard
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

global Pq cq Vq
rhsp = Pq*(cq.*(Vq*rhsp));


return;

