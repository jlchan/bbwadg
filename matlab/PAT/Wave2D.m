function Wave2D

% clear all, clear
clear -global *

Globals2D

N = 5;
K1D = 32;
c_flag = 0;
cfun = @(x,y) ones(size(x));
%     cfun = @(x,y) 1 + .5*sin(pi*x).*sin(pi*y); % smooth velocity
% cfun = @(x,y) (1 + .5*sin(2*pi*x).*sin(2*pi*y) + (y > 0)); % piecewise smooth velocity

% filename = 'Grid/Other/block2.neu';
% filename = 'Grid/Maxwell2D/Maxwell05.neu';
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
StartUp2D;

BuildPeriodicMaps2D(1,1);

global Pq cq Vq

% PlotMesh2D; return

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

%%
cq = cfun(xq,yq);

%FinalTime = 2/min(cq);
FinalTime = 1/min(cq);

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

% p = cos(W*x).*cos(W*y);
x0 = 0; y0 = .1;
% p = exp(-200*((x-x0).^2 + (y-y0).^2));
p = -exp(-30^2*((x-x0).^2 + (y-y0).^2)); 
u = zeros(Np, K); 
v = zeros(Np, K); 

%%

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
CNh = max(CN*max(Fscale(:)));
dt = 2/CNh;

% outer time step loop
tstep = 0;
figure

while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
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
    
    if 1 && nargin==0 && mod(tstep,10)==0
        clf
        pp = p;
        vv = Vp*pp;        
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
%         caxis([-.1 .2])
%         axis([0 1 0 1 -10 10])
        %         PlotField2D(N+1, x, y, pp); view(2)
        title(sprintf('time = %f',time))
        drawnow
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
end

axis off
view(0,0)


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

% % Impose reflective boundary conditions (p+ = -p-)
% ndotdU(mapB) = 0;
% dp(mapB) = -2*p(vmapB);

% % free surface BCs
% ndotdU(mapB) = -2*(nx(mapB).*u(vmapM(mapB)) + ny(mapB).*v(vmapM(mapB)));
% dp(mapB) = 0;

% basic ABCs
fluxp(mapB) = (nx(mapB).*u(vmapM(mapB)) + ny(mapB).*v(vmapM(mapB)));
fluxu(mapB) = p(vmapM(mapB)).*nx(mapB);
fluxv(mapB) = p(vmapM(mapB)).*ny(mapB);

tau = 1;
fluxp =  tau*dp - ndotdU;
fluxu =  (tau*ndotdU - dp).*nx;
fluxv =  (tau*ndotdU - dp).*ny;

pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
rhsp =  -divU + LIFT*(Fscale.*fluxp)/2.0;
rhsu =  -dpdx + LIFT*(Fscale.*fluxu)/2.0;
rhsv =  -dpdy + LIFT*(Fscale.*fluxv)/2.0;

global Pq cq Vq
rhsp = Pq*(cq.*(Vq*rhsp));



return;

