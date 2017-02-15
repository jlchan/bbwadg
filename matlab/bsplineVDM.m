function [BVDM M Dr] = bsplineVDM(N,K,rin)

if nargin==0
    N = 4;
    K = 4;
    rin = JacobiGL(0,0,N+K-1);
end

[Nv, VX, K, EToV] = MeshGen1D(-1,1,K);

map = @(r) reshape((1/K)*(repmat(r,1,K)+1) + ...
    repmat(VX(1:end-1),length(r),1),length(r)*K,1);

% extraction nodes - interpolation
r = JacobiGL(0,0,N);
rB = map(r);

% local quadrature
[rq, wq] = JacobiGQ(0,0,2*N+1);
rBq = map(rq);
wBq = repmat(wq,1,K)*(1/K);
wBq = wBq(:);

% extract local basis
t = [VX(1)*ones(1,N) VX VX(end)*ones(1,N)]; % open knot vec
R = bspline_basismatrix(N+1,t,rB); % for local extraction

% local basis
V = Vandermonde1D(N,r);
Drq = GradVandermonde1D(N,rq)/V;

Bq = bspline_basismatrix(N+1,t,rBq);
M = Bq'*diag(wBq)*Bq;
S = Bq'*diag(wBq)*kron(eye(K)*K,Drq)*R; % rx = K1D in 1D
Dr = M\S;

% VDM for interp or whatever
BVDM = bspline_basismatrix(N+1,t,rin);

% rp = linspace(-1,1,250);
% Vp = Vandermonde1D(N,rp)/V;
% Bp = bspline_basismatrix(N+1,t,rp);
% plot(rp,Bp,'.--')
