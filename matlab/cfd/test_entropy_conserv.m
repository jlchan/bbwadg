clear
Globals1D;

N = 7;
K1D = 2;

[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
StartUp1D;
vmapP(1) = vmapM(end); vmapP(end) = vmapM(1);
mapM = reshape(1:Nfp*Nfaces*K,Nfp*Nfaces,K);
mapP = mapM;
for e = 1:K
    mapP(:,e) = [(Nfp*Nfaces)*e-2, Nfp*Nfaces*(e+1)-1];
end
mapM = mapM(:); mapP = mapP(:);
mapP(1) = mapM(end); mapP(end) = mapM(1);

r = JacobiGL(0,0,N);
% r = JacobiGQ(0,0,N);
[rq wq] = JacobiGQ(0,0,N);
[rq wq] = JacobiGQ(0,0,N+1);

% % % include boundary nodes
% rq = [-1;rq;1];
% wq = [0;wq;0];
V = Vandermonde1D(N,r);

Nq = length(rq);
Vq = Vandermonde1D(N,rq)/V;

% make global points
rp = linspace(-1,1,100)';
Vp = Vandermonde1D(N,rp)/V;
xq =  Vandermonde1D(N,rq)/Vandermonde1D(N,JacobiGL(0,0,N))*x;
xp = Vandermonde1D(N,rp)/Vandermonde1D(N,JacobiGL(0,0,N))*x;
x = x(:); xq = xq(:); xp = xp(:);
wJq = diag(wq)*(Vq*J); wJq = wJq(:);

% usual local ops
V = Vandermonde1D(N,r);
Dr = GradVandermonde1D(N,r)/V;
Vf = Vandermonde1D(N,[-1 1])/V;
M = Vq'*diag(wq)*Vq;
LIFT = M\Vf';
Pq = M\(Vq'*diag(wq));
Prq = M\((Vq*Dr)'*diag(wq));
W = diag(wq);

Sq = W*Vq*Dr*Pq;
B = Pq'*Vf'*diag([-1,1])*Vf*Pq;

Sw = B-Pq'*Dr'*Vq'*diag(wq);

uq = exp(cos(1+2*rq));
% uq = exp(cos(1+rq));

[ux uy] = meshgrid(uq);
F = logmean(ux,uy); v = log(uq);

% Dskew = Vq*Dr*Pq - Vq*(M\(Dr'*Vq'*W));

r1 = v'*sum(2*Sq.*F,2);
uf = Vf*Pq*uq;
r2 = v'*sum(B.*F,2) - (uf(2)-uf(1));
r3 = (Vq*Pq*v)'*sum((2*Sq).*F,2); % testing with projection of V
% r4 = v'*W*Vq*Pq*sum(2*Dq.*F,2); % projecting flux before testing
norm(r1-r2)
norm(r2-r3)

%(Vq*Pq*v)'*sum(B.*F,2) - (uf(2)-uf(1))
% ux-uy
% F = (ux+uy)/2; v = uq;

% sum(S*uq)
% sum(B*uq)


