
Globals1D;

NB = 4; Ksub = 8; K1D = 1;
% NB = 4; Ksub = 1; K1D = 4;
% NB = 11; Ksub = 1; K1D = 1;
smoothKnots = 25;
N = NB+Ksub-1;

FinalTime = .5;

[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

rp = linspace(-1,1,250);
Vp = Vandermonde1D(N,rp)/V;
xp = Vp*x;
xe = Vandermonde1D(N,linspace(-1,1,N+1)')/V*x;
dofs = (N+1)*K1D;

% make splines
rB = JacobiGL(0,0,N);
xB = Vandermonde1D(N,rB)/V * x;

[BVDM M Dr R rq wq Vq, ~, rx] = bsplineVDM(NB,Ksub,rB,smoothKnots,N); % VDM for interp, mass, M\S

xq = (Vandermonde1D(N,rq)/V)*x;
wJq = diag(wq)*(Vq*J);

Vf = zeros(2,N+1);
Vf(1,1) = 1;
Vf(2,N+1) = 1;
% keyboard
Lq = M\Vf';

[Vp] = bsplineVDM(NB,Ksub,rp,smoothKnots);

Pq = M\(Vq'*diag(wq));

f = @(x) exp(-100*(x-.2).^2);
f = @(x) x > .3;
plot(xp,f(xp),'linewidth',2)
hold on
plot(xp,Vp*Pq*f(xq),'--','linewidth',2)
% plot(xe,Pq*f(xq),'o--','linewidth',2)
% plot(xp,Vp*f(xe),'--','linewidth',2)

VqPq = Vq*Pq;
Drq = Vq*Dr*Pq;

VqLq = Vq*Lq;
VfPq = Vf*Pq;

nrJ = [-1;1];
DNr = [Drq-.5*Vq*Lq*diag(nrJ)*VfPq .5*VqLq*diag(nrJ); -.5*diag(nrJ)*VfPq .5*diag(nrJ)*eye(2)];
W = diag([wq;1;1]);


