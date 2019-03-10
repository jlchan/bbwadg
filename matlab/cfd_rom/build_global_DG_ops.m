clear
clear -global
Globals1D;

N = 1;
K1D = 10;

[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

global QNr Nq

[r w] = JacobiGL(0,0,N);
% [rq wq] = JacobiGQ(0,0,N+1);
[rq wq] = JacobiGQ(0,0,N);
Nq = length(rq);

V = Vandermonde1D(N,r);
Dr = GradVandermonde1D(N,r)/V;
Vq = Vandermonde1D(N,rq)/V;
Vf = Vandermonde1D(N,[-1 1])/V;

M = Vq'*diag(wq)*Vq;
Lq = M\Vf';
Pq = M\(Vq'*diag(wq));
Pq(abs(Pq)<1e-8) = 0;

VqPq = Vq*Pq;
VqLq = Vq*Lq;
VfPq = Vf*Pq;
Drq = Vq*Dr*Pq;
nrJ = [-1;1];
Qr = diag(wq)*Vq*Dr*Pq;
Qrskew = (Qr-Qr');
QNr = .5*[Qrskew (Vf*Pq)'*diag([-1;1])
    -diag([-1;1])*Vf*Pq zeros(2)];

E = (Vf*Pq);
B = diag([-1;1]);
QNr = [Qr-.5*E'*B*E .5*E'*B
    -.5*B*E .5*B];

WN = diag([wq;1;1]);
PN = M\[Vq' Vf']*WN;
DNr = WN\QNr;
W = diag([wq;1;1]);
Ef = zeros(2,length(rq)); Ef(1) = 1; Ef(end) = 1;

wJq = diag(wq)*(Vq*J);

% maps
mapM = reshape(1:Nfp*Nfaces*K,Nfp*Nfaces,K);
mapP = mapM;
for e = 1:K
    mapP(:,e) = [(Nfp*Nfaces)*e-2, Nfp*Nfaces*(e+1)-1];
end
mapP(1) = 1; mapP(end) = Nfp*Nfaces*K;
mapP = reshape(mapP,Nfp*Nfaces,K);
mapP(1) = mapM(end); mapP(end) = mapM(1); % periodic

u = zeros(Nq+2,K);
Q = zeros((Nq+2)*K);
Bx = zeros(Nq+2,2);
Bx(Nq+1:Nq+2,:) = diag([-1;1]);

for i = 1:(Nq+2)*K     
    u(i) = 1;    
    uf = u(Nq+1:Nq+2,:);
    rhs = QNr*u + .5*Bx*(uf(mapP)-uf); 
    u(i) = 0;
    
    Q(:,i) = rhs(:);           
end
ids = reshape(1:(Nq+2)*K,Nq+2,K);
idq = ids(1:Nq,:);
idf = ids(Nq+1:Nq+2,:);
p = [idq(:); idf(:)];
Q = Q(p,p);

idq = 1:Nq*K;
idf = Nq*K+1:(Nq+2)*K;
% spy(Q)

VqK = kron(eye(K),Vq);
VfK = kron(eye(K),Vf);

MK = VqK'*diag(wJq(:))*VqK;
PqK = MK\(VqK'*diag(wJq(:)));
PK = MK\[VqK;VfK]';

% reduced basis
xq = Vq*x;
xf = Vf*x;
xN = [VqK;VfK]*x(:);

Vr = PqK*orth([xq(:).^(0:2)]);
% Vr = PqK*orth([xq(:).^0 sin(pi*xq(:)) cos(pi*xq(:)) sin(2*pi*xq(:)) cos(2*pi*xq(:)) sin(3*pi*xq(:)) cos(3*pi*xq(:)) sin(4*pi*xq(:)) cos(4*pi*xq(:))]);
% Vr = eye((N+1)*K);

VB = VfK([1 end],:);
f = @(x) exp(sin(pi*x));
df = @(x) pi*cos(pi*x).*exp(sin(pi*x));

Pr = (Vr'*MK*Vr)\(Vr'*VqK'*diag(wJq(:)));
Er = VB*Vr*Pr;
QB = [(kron(eye(K),Qr)*VqK*Vr*Pr-.5*Er'*B*Er) .5*Er'*B; -.5*B*Er .5*B];
plot(x(:),Vr*((Vr'*MK*Vr)\(Vr'*[VqK;VfK]'*Q*f(xN(:))))-df(x(:)),'bo--','linewidth',2,'markersize',12)
hold on
plot(x(:),Vr*((Vr'*MK*Vr)\(Vr'*[VqK;VB]'*QB*f([xq(:);-1;1])))-df(x(:)),'rx--','linewidth',2,'markersize',12)
% plot(x(:),Vr*((Vr'*MK*Vr)\(Vr'*VqK'*kron(eye(K),Qr)*f(xq(:))))-df(x(:)),'g^--','linewidth',2,'markersize',12)
% plot(x(:),df(x(:)))
PK = MK\[VqK;VfK]';

% plot(x(:),PK*Q*f(xN),'o--')
% hold on
% plot(x(:),f(x(:)),'-')


