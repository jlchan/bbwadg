function Wave2D

% clear all, clear
clear -global *

Globals2D

N = 7;
K1D = 8;
c_flag = 0;
cfun = @(x,y) ones(size(x));
% cfun = @(x,y) 1 + (x > 0);
%     cfun = @(x,y) 1 + .5*sin(pi*x).*sin(pi*y); % smooth velocity
% cfun = @(x,y) (1 + .5*sin(2*pi*x).*sin(2*pi*y) + (y > 0)); % piecewise smooth velocity

% filename = 'Grid/Other/block2.neu';
% filename = 'Grid/Maxwell2D/Maxwell05.neu';
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
% VX = 1*(K1D)/(K1D-2)*VX; VY = 1*(K1D)/(K1D-2)*VY;
% VX = 2*VX; VY = 2*VY;
StartUp2D;

% BuildPeriodicMaps2D(1,1);

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

% FinalTime = 2/min(cq(:));
FinalTime = (max(VX)-min(VX))/min(cq(:))

FinalTime = .3;
%%

global mapBx vmapBx
%
% % find x = 0 faces
% fbids = reshape(find(vmapP==vmapM),Nfp,nnz(vmapP==vmapM)/Nfp);
% xfb = x(vmapM(fbids));
% bfaces = sum(abs(xfb+1),1)<1e-8;
%
% fbids = fbids(:,bfaces);
% fbids = fbids(:);
%
% mapBx = fbids;
% vmapBx = vmapM(mapBx);

% get external traces
vmapMx = reshape(vmapM,Nfp*Nfaces,K);
% vmapMx = vmapMx(:,sum(x < -1,1)>0);
vmapBx = vmapMx(abs(x(vmapMx)+1)<1e-8 & abs(y(vmapMx))<1+1e-8 );
% plot(x(vmapMx),y(vmapMx),'s','markersize',10);
% return

%% params setup

% p = cos(W*x).*cos(W*y);
x0 = .1; y0 = 0;
% p = exp(-200*((x-x0).^2 + (y-y0).^2));
p = exp(-50^2*((x-x0).^2 + (y-y0).^2));
p = Pq*(abs(xq)<.5 & abs(yq)<.5);

sig = @(x) 1-1./(1+exp(-250*x));
r2 = @(x,y) x.^2 + y.^2;
p = sig((r2(x,y)-.25).^2);

u = zeros(Np, K);
v = zeros(Np, K);

% vv = Vp*p; color_line3(xp,yp,vv,vv,'.'); return

%% check eigs

if 3*Np*K < 4000
    U = zeros(Np*K,3);
    A = zeros(3*Np*K);
    for i = 1:3*Np*K
        U(i) = 1;
        p = reshape(U(:,1),Np,K);
        u = reshape(U(:,2),Np,K);
        v = reshape(U(:,3),Np,K);
        [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,0);
        A(:,i) = [rhsp(:); rhsu(:); rhsv(:)];
        U(i) = 0;
        if mod(i,round(3*Np*K/10))==0
            disp(sprintf('i = %d out of %d\n',i,3*Np*K))
        end
    end
    lam = eig(A);
    max(abs(lam))
    keyboard
    
end

%%

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
CNh = max(CN*max(Fscale(:)));
dt = 3/CNh
Nstep = ceil(FinalTime/dt);
dt = FinalTime/Nstep
phist = zeros(length(vmapBx),ceil(FinalTime/dt));

% outer time step loop
tstep = 0;
figure
time = 0;
for i = 1:Nstep
    
    for INTRK = 1:5
        
        timelocal = time + rk4c(INTRK)*dt;
        
        [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,timelocal);
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        p = p+rk4b(INTRK)*resp;
        
    end;
    
    if 1 && nargin==0 && mod(tstep,10.03512)==0
        clf
        vv = Vp*p;
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        
%         ids = abs(yp)<1e-8;
%         plot(xp(ids),vv(ids),'.')
%         axis([-1 1 -.15 .15])

%         caxis([-.25 .26])
        %         axis([0 1 0 1 -10 10])
        %         PlotField2D(N+1, x, y, pp); view(2)
        title(sprintf('time = %f, max solution val = %f',time,max(abs(vv(:)))))
        drawnow
    end
    
    % Increment time
    time = time+dt; 
    
    phist(:,i) = p(vmapBx);
    thist(:,i) = time;
end

keyboard

%%
figure
plot3(repmat(y(vmapBx),1,tstep),thist,phist,'.'); xlabel('x')



function [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

% evaluate upwind fluxes
ndotdU = nx.*du + ny.*dv;

% Impose reflective boundary conditions (p+ = -p-)
ndotdU(mapB) = 0;
dp(mapB) = -2*p(vmapB);

% % free surface BCs
% ndotdU(mapB) = -2*(nx(mapB).*u(vmapM(mapB)) + ny(mapB).*v(vmapM(mapB)));
% dp(mapB) = 0;

% % % basic abcs
% dp(mapB) = -p(vmapB);
% ndotdU(mapB) = -(nx(mapB).*u(vmapM(mapB)) + ny(mapB).*v(vmapM(mapB)));

tau = 1;
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

