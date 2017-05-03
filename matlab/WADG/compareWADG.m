function Wave2D

% clear all, clear
clear -global *

Globals2D

N = 5;
K1D = 8;
c_flag = 0;
FinalTime = 50;
cfun = @(x,y) ones(size(x));
cfun = @(x,y) 1 + .5*sin(pi*x).*sin(pi*y); % smooth velocity
% cfun = @(x,y) (1 + .5*sin(2*pi*x).*sin(2*pi*y) + (y > 0)); % piecewise smooth velocity

filename = 'Grid/Other/block2.neu';
% filename = 'Grid/Maxwell2D/Maxwell05.neu';
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
% VX = (VX + 1)/2;
% VY = (VY + 1)/2;
StartUp2D;

% BuildPeriodicMaps2D(1,1);

% PlotMesh2D; return
%%

global Vq Pq cq wq

[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
Pq = V*V'*Vq'*diag(wq); % J's cancel out
Mref = inv(V*V');
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;
wJq = diag(wq)*Jq;

cq = cfun(xq,yq);

% color_line3(xq,yq,cq,cq,'.');

%%

global mapBx vmapBx

% find x = 0 faces
fbids = reshape(find(vmapP==vmapM),Nfp,nnz(vmapP==vmapM)/Nfp);
xfb = x(vmapP(fbids));
bfaces = sum(abs(xfb),1)<1e-8;

fbids = fbids(:,bfaces);
fbids = fbids(:);

mapBx = fbids;
vmapBx = vmapP(mapBx);


%% params setup

global t0
t0 = .1;

k = 1; % frequency of solution
W = (2*k-1)/2*pi;
% p = cos(W*x).*cos(W*y);
x0 = 0; y0 = 0;
p = exp(-50*((x-x0).^2 + (y-y0).^2));
% p = zeros(Np, K); 
u = zeros(Np, K); 
v = zeros(Np, K);

p2 = p;
u2 = u;
v2 = v;

%%

if 0
    A = zeros(3*Np*K);
    U = zeros(3*Np,K);
    for i = 1:3*Np*K
        U(i) = 1;
        ids = 1:Np;
        p = U(ids,:);
        u = U(Np+ids,:);
        v = U(2*Np+ids,:);
        [rhsp, rhsu, rhsv] = acousticsRHS2D_DG(p,u,v,0);
        rhs = [rhsp; rhsu; rhsv];
        A(:,i) = rhs(:);
        U(i) = 0;
        if mod(i,round(3*Np*K/10))==0
            disp(sprintf('i = %d out of %d\n',i,3*Np*K));
        end
    end
    lam = eig(A);
    plot(lam,'o')
    keyboard
end
%%

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);
resu2 = zeros(Np,K); resv2 = zeros(Np,K); resp2 = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
dt = 1/max(CN*max(Fscale(:)));

% outer time step loop
tstep = 0;
figure
colormap(gray)
while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        timelocal = time + rk4c(INTRK)*dt;
        
        [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,timelocal);
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        p = p+rk4b(INTRK)*resp;
        
        [rhsp, rhsu, rhsv] = acousticsRHS2D_DG(p,u,v,timelocal);
        resp2 = rk4a(INTRK)*resp2 + dt*rhsp;
        resu2 = rk4a(INTRK)*resu2 + dt*rhsu;
        resv2 = rk4a(INTRK)*resv2 + dt*rhsv;
        u2 = u2+rk4b(INTRK)*resu2;
        v2 = v2+rk4b(INTRK)*resv2;
        p2 = p2+rk4b(INTRK)*resp2;
        
    end;
    
    if 1 && nargin==0 && mod(tstep,10)==0
        clf
        subplot(1,2,1)
        vv = Vp*p;
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        
        subplot(1,2,2)
        vv = Vp*p2;
        color_line3(xp,yp,vv,vv,'.');        
        axis equal
        axis tight
        
%         colorbar
%         caxis([-.1 .2])
%         axis([0 1 0 1 -10 10])
        %         PlotField2D(N+1, x, y, pp); view(2)
        title(sprintf('time = %f',time))
        drawnow
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
    
    tvec(tstep) = time;
    udiff = (Vq*(p-p2)).^2;
    diff(tstep) = sqrt(sum(wJq(:).*udiff(:)));
end

figure
semilogy(tvec,diff,'.-','linewidth',2)

keyboard

% axis off
% view(0,0)


function [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

global Pq Vq cq

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

% evaluate upwind fluxes
ndotdU = nx.*du + ny.*dv;

% % % Impose reflective boundary conditions (p+ = -p-)
% ndotdU(mapB) = 0;
% dp(mapB) = -2*p(vmapB);

% free surface BCs
ndotdU(mapB) = -2*(nx(mapB).*u(vmapM(mapB)) + ny(mapB).*v(vmapM(mapB)));
dp(mapB) = 0;

tau = 1;
fluxp =  tau*dp - ndotdU;
fluxu =  (tau*ndotdU - dp).*nx;
fluxv =  (tau*ndotdU - dp).*ny;

% % basic ABCs
% fluxp(mapB) = (nx(mapB).*u(vmapM(mapB)) + ny(mapB).*v(vmapM(mapB)));
% fluxu(mapB) = p(vmapM(mapB)).*nx(mapB);
% fluxv(mapB) = p(vmapM(mapB)).*ny(mapB);

pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
rhsp =  -divU + LIFT*(Fscale.*fluxp)/2.0;
rhsu =  -dpdx + LIFT*(Fscale.*fluxu)/2.0;
rhsv =  -dpdy + LIFT*(Fscale.*fluxv)/2.0;

rhsp = Pq*(cq.*(Vq*rhsp));


function [rhsp, rhsu, rhsv] = acousticsRHS2D_DG(p,u,v,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

global Vq Pq cq wq

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

% evaluate upwind fluxes
ndotdU = nx.*du + ny.*dv;

% % % Impose reflective boundary conditions (p+ = -p-)
% ndotdU(mapB) = 0;
% dp(mapB) = -2*p(vmapB);

% free surface BCs
ndotdU(mapB) = -2*(nx(mapB).*u(vmapM(mapB)) + ny(mapB).*v(vmapM(mapB)));
dp(mapB) = 0;

tau = 1;
fluxp =  tau*dp - ndotdU;
fluxu =  (tau*ndotdU - dp).*nx;
fluxv =  (tau*ndotdU - dp).*ny;

% % basic ABCs
% fluxp(mapB) = (nx(mapB).*u(vmapM(mapB)) + ny(mapB).*v(vmapM(mapB)));
% fluxu(mapB) = p(vmapM(mapB)).*nx(mapB);
% fluxv(mapB) = p(vmapM(mapB)).*ny(mapB);

pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
rhsp =  -divU + LIFT*(Fscale.*fluxp)/2.0;
rhsu =  -dpdx + LIFT*(Fscale.*fluxu)/2.0;
rhsv =  -dpdy + LIFT*(Fscale.*fluxv)/2.0;

for e = 1:K
    Mc = V*V'*Vq'*diag(wq(:)./cq(:,e))*Vq;
    rhsp(:,e) = Mc\rhsp(:,e);
end


