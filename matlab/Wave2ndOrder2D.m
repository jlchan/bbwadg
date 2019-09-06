clear all -globals
Globals2D

N = 5;
K1D = 8;
FinalTime = .50;
CFL = .9;

global tau
tau = 1;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
StartUp2D;
BuildPeriodicMaps2D(max(VX)-min(VX),max(VY)-min(VY));

vmapM = reshape(vmapM,Nfp*Nfaces,K);
vmapP = reshape(vmapP,Nfp*Nfaces,K);
mapP  = reshape(mapP, Nfp*Nfaces,K);
nx = reshape(nx,Nfp*Nfaces,K);
ny = reshape(ny,Nfp*Nfaces,K);
sJ = reshape(sJ,Nfp*Nfaces,K);

% plotting nodes
[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

global Pq Vq
Nq = 2*N;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
M = Vq'*diag(wq)*Vq; % I prefer this approach since M != inv(V*V) if quadrature is inexact
Pq = M\(Vq'*diag(wq)); % J's cancel out
xq = Vq*x; yq = Vq*y;

global rxJ sxJ ryJ syJ nrJ nsJ

rxJ = rx.*J; ryJ = ry.*J;
sxJ = sx.*J; syJ = sy.*J;

nrJ = zeros(Nfp,Nfaces);
nsJ = zeros(Nfp,Nfaces);
nrJ(:,1) = 0; nrJ(:,2) = 1; nrJ(:,3) = -1;
nsJ(:,1) = -1; nsJ(:,2) = 1; nsJ(:,3) = 0;

nrJ = nrJ(:);
nsJ = nsJ(:);

%% params setup

x0 = 0;
y0 = 0;
p = exp(-10^2*((x-x0).^2 + (y-y0).^2));

global cq
cq = 1 + .25*sin(pi*xq).*sin(pi*yq);

%% compute eigs

if 0
    p = zeros(Np,K);
    for i = 1:Np*K
        p(i) = 1;
        rhs = rhs2ndOrder(p);
        A(:,i) = rhs(:);
        p(i) = 0;
    end
    M = kron(diag(J(1,:)),Vq'*diag(wq)*Vq);
    keyboard
end

%%

pprev = p;

CN = (N+1)*(N+2)/2; % guessing...
CNh = max(CN*max(Fscale(:)));
dt = CFL*2/CNh;

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

for i = 1:Nsteps
    
    % leapfrog: p(k+1) = 2*p(k) - p(k-1) + dt^2*rhs = 0;
    [p,pprev] = rhs2ndOrder2(p,pprev,dt);   
        
    if mod(i,10)==0 || i==Nsteps
        clf
        vv = Vp*p;
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('step = %d/%d, time = %f',i,Nsteps,i*dt))
        drawnow
    end
        
end

% simple straightforward version with physical geofacs
function [p pprev] = rhs2ndOrder(p,pprev,dt)

Globals2D;

% compute DG gradient
dp = p(vmapP)-p(vmapM);
pr = Dr*p; 
ps = Ds*p;
px = rx.*pr + sx.*ps + .5*LIFT*(Fscale.*dp.*nx);
py = ry.*pr + sy.*ps + .5*LIFT*(Fscale.*dp.*ny);

dpx = px(vmapP)-px(vmapM);
dpy = py(vmapP)-py(vmapM);
dpn = dpx.*nx + dpy.*ny;

% compute DG divergence
pxr = Dr*px; pxs = Ds*px;
pyr = Dr*py; pys = Ds*py;
global tau
rhs = rx.*pxr + sx.*pxs + ry.*pyr + sy.*pys + .5*LIFT*(Fscale.*(dpn-tau*dp));

% update
pnew = 2*p - pprev + dt^2*rhs;
pprev = p;
p = pnew;

end

% a refactored version using only volume geofacs
function [p pprev] = rhs2ndOrder2(p,pprev,dt)

Globals2D;
global rxJ sxJ ryJ syJ nrJ nsJ

% compute DG gradient
dp = p(vmapP)-p(vmapM);
pr = Dr*p + .5*LIFT*(diag(nrJ)*dp); 
ps = Ds*p + .5*LIFT*(diag(nsJ)*dp);

% rotate to physical element
px = rx.*pr + sx.*ps;
py = ry.*pr + sy.*ps;

% not sure if I can move this part, since geofacs are discontinuous
dpx = px(vmapP)-px(vmapM);
dpy = py(vmapP)-py(vmapM);

% compute DG divergence 
pr = rx.*px + ry.*py; 
ps = sx.*px + sy.*py; 
dpr = rx(Fmask(:),:).*dpx + ry(Fmask(:),:).*dpy;
dps = sx(Fmask(:),:).*dpx + sy(Fmask(:),:).*dpy;
dpn = diag(nrJ)*dpr + diag(nsJ)*dps;

global tau
rhs = Dr*pr + Ds*ps + .5*LIFT*(dpn - tau*dp);

% update
pnew = 2*p - pprev + dt^2*rhs;
pprev = p;
p = pnew;

end
