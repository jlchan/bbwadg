function Wave2D

clear all -globals
Globals2D

N = 7;
K1D = 8;
c_flag = 0;
FinalTime = .50;
CFL = .5;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
StartUp2D;

% plotting nodes
[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

global Pq Vq
Nq = 2*N;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
M = Vq'*diag(wq)*Vq; % I prefer this approach since M != inv(V*V) if quadrature is inexact
Pq = M\(Vq'*diag(wq)); % J's cancel out
Mref = inv(V*V');
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

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
t0 = .0;

k = 1; % frequency of solution
W = (2*k-1)/2*pi;
p = cos(W*x).*cos(W*y);
x0 = 0; y0 = .1;
p = exp(-200*((x-x0).^2 + (y-y0).^2));

x0 = 0;
y0 = .25;
p = exp(-10^2*((x-x0).^2 + (y-y0).^2));
% p = zeros(Np, K); 
u = zeros(Np, K); 
v = zeros(Np, K);

% v = -exp(-30^2*((x-x0).^2 + (y-y0).^2));

% x0 = 0;
% y0 = -.25;
% v = -exp(-10^2*((x-x0).^2 + (y-y0).^2));
% p = exp(-50^2*((x-x0).^2));
% u = zeros(Np, K);


%%
x0 = -.1; y0 = 0;
for e = 1:K
    cx(e) = mean(VX(EToV(e,:)));
    cy(e) = mean(VY(EToV(e,:)));    
end

[~,emin] = min(abs(cx - x0) + abs(cy-y0));
x0 = cx(emin); y0 = cy(emin);

global rick ptsrc 
f0 = 170;
f0 = 100;
tR = 1/f0;
rick = @(t) 1e6*(1 - 2*(pi*f0*(t-tR)).^2).*exp(-(pi*f0*(t-tR)).^2).*(t < tR);
% tt = 0:.0001:.1; plot(tt,rick(tt));return
ptsrc = @(x,y) exp(-(200)^2*((x-x0).^2 + (y-y0).^2));

if 1
    [rq2 sq2 wq2] = Cubature2D(4*N+2);
    Vq2 = Vandermonde2D(N,rq2,sq2)/V;
    Pq2 = (V*V') * Vq2'*diag(wq2);
    xq2 = Vq2*x; yq2 = Vq2*y;
    ptsrc = ptsrc(xq2,yq2);
    err = max(max(abs(ptsrc - Vq2*Pq2*ptsrc)))
    ptsrc = Pq2*ptsrc;
    ptsrc = ptsrc/max(abs(ptsrc(:)));
else % interp
    err = max(max(abs(Vp*ptsrc(x,y) - ptsrc(xp,yp))))
    ptsrc = ptsrc(x,y);    
end

ptsrc(:) = 0; 
% ptsrc(:,emin) = (V*V')*(Vandermonde2D(N,mean([-1 -1 1]), mean([-1 -1 1]))/V)'; % ptwise evaluation

%%

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
CNh = max(CN*max(Fscale(:)));
dt = CFL*2/CNh;

% outer time step loop
tstep = 0;
figure
% colormap(gray)
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
        %         vv = abs(vv);
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

keyboard
% axis off
% view(0,0)


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

% left wall bcs
global mapBx vmapBx t0
if (time < t0)
    dp(mapBx) = 2*(sin(pi*time/t0)*(time<t0) - p(vmapBx));
    ndotdU(mapBx(:)) = 0;
end

% % basic ABCs
% fluxp(mapB) = (nx(mapB).*u(vmapM(mapB)) + ny(mapB).*v(vmapM(mapB)));
% fluxu(mapB) = p(vmapM(mapB)).*nx(mapB);
% fluxv(mapB) = p(vmapM(mapB)).*ny(mapB);

tau = 1;
fluxp =  tau*dp - ndotdU;
fluxu =  (tau*ndotdU - dp).*nx;
fluxv =  (tau*ndotdU - dp).*ny;


pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
global rick ptsrc 
rhsp =  -divU + LIFT*(Fscale.*fluxp)/2.0 + rick(time)*ptsrc;
rhsu =  -dpdx + LIFT*(Fscale.*fluxu)/2.0;
rhsv =  -dpdy + LIFT*(Fscale.*fluxv)/2.0;

global Pq cq Vq
% rhsp = Pq*(cq.*(Vq*rhsp));
rhsp = rhsp.*(Pq*cq);



return;

