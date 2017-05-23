function [BVDM M Dr R rBq wBq Bq Brq VX] = bsplineVDM(N,K,rin,smoothKnots)
 
if nargin==0
    N = 4;
    K = 4;
    rin = JacobiGL(0,0,N+K-1);    
end

if nargin < 4
    smoothKnots = 0;
end

if (K==1)
    BVDM = bern_basis_1D(N,rin);
    [rBq wBq] = JacobiGQ(0,0,N);
    [Vq Vrq] = bern_basis_1D(N,rBq);
    M = Vq'*diag(wBq)*Vq;
    Dr = M\(Vq'*diag(wBq)*Vrq);
    R = eye(N+1);
    Bq = Vq;
    Brq = Vrq;
    VX = [-1 1];
    return
end

[Nv, VX, K, EToV] = MeshGen1D(-1,1,K);

map = @(r) reshape((1/K)*(repmat(r,1,K)+1) + ...
    repmat(VX(1:end-1),length(r),1),length(r)*K,1);

t0 = [VX(1)*ones(1,N) VX VX(end)*ones(1,N)];
t = t0;
if strcmp(smoothKnots,'opt')
    
    if N==1
        VX = linspace(-1,1,Ksub+1); % optimal for N=1
    elseif K==2*N || K==N
        load optKnots.mat        
        VX = optKnots{NBlist==N & KsubList==K};        
    else        
        if N==4 && K==16
            VX = [-1.0000   -0.8004   -0.6856   -0.5714   -0.4571   -0.3429   -0.2286   -0.1143         0    0.1143    0.2286    0.3429    0.4571    0.5714    0.6856    0.8004    1.0000];
        elseif N==4 && K==32
            VX = [-1.0000   -0.8957   -0.8357   -0.7761   -0.7164   -0.6567   -0.5970   -0.5373   -0.4776   -0.4179   -0.3582   -0.2985   -0.2388   -0.1791   -0.1194   -0.0597         0    0.0597    0.1194    0.1791    0.2388    0.2985   0.3582    0.4179    0.4776    0.5373    0.5970    0.6567    0.7164    0.7761    0.8357    0.8957    1.0000];
        end

    end
    t0 = [VX(1)*ones(1,N) VX VX(end)*ones(1,N)];
    t = t0;
    
elseif smoothKnots > 0
    re = linspace(-1,1,N+K)';
    for i = 1:smoothKnots        
        Ve = bspline_basismatrix(N+1, t, t0);
        t = Ve*re;        
    end
    VX = t(N+1:end-N)';
    
end

h = @(r) repmat(diff(VX),length(r),1);
map = @(r) reshape(h(r).*(repmat(r,1,K)+1)/2 + repmat(VX(1:end-1),length(r),1),length(r)*K,1);

% re = linspace(-1,1,N+1)';
% Ve = bspline_basismatrix(N+1,[VX(1)*ones(1,N) VX VX(end)*ones(1,N)],re);
% VX = Ve*re;

% extraction nodes - interpolation
r = JacobiGL(0,0,N);
rB = map(r);

% local quadrature
if 1
    [rq, wq] = JacobiGQ(0,0,N);
    rBq = map(rq);
    rBq = rBq(:);
    wBq = repmat(wq,1,K).*h(rq)/2;    
    wBq = wBq(:);
else
    [rq, wq] = JacobiGQ(0,0,N-2);
    rBq = map(rq);
    rBq = rBq(:);
    wBq = repmat(wq,1,K).*h(rq)/2;
    wBq = wBq(:);
    [rBq,ids] = uniquetol(rBq);
    wBq = wBq(ids);
%     keyboard    
end


% extract local basis
% t = [VX(1)*ones(1,N) VX VX(end)*ones(1,N)]; % open knot vec
R = bspline_basismatrix(N+1,t,rB); % for local extraction to Lagrange dofs
R(abs(R)<1e-8) = 0;
Bq = bspline_basismatrix(N+1,t,rBq);
M = Bq'*diag(wBq)*Bq;
    
% local basis
V = Vandermonde1D(N,r);

if 1
    % local projection
    Drq = GradVandermonde1D(N,rq)/V;    
    Brq = kron(diag(2./diff(VX)),Drq)*R;    
    
    S = Bq'*diag(wBq)*Brq; % rx = 2/h in 1D
    Dr = M\S;    
    
else % non-local quadratures
    
    % spline derivative matrices
    dt = [VX(1)*ones(1,N) sort([VX(2:end-1) VX(2:end-1)]) VX(end)*ones(1,N)];
    DB = bspline_basismatrix(N,dt,rB); % for local extraction to Lagrange dofs

    % local derivatives
    Dr = GradVandermonde1D(N,r)/V;
    Br = kron(diag(2./diff(VX)),Dr)*R; 
    
    % lsq fit 
    DBr = DB\Br;
    DBr = sparse(DBr);
    keyboard
    
    DBq = bspline_basismatrix(N,dt,rBq);
        
    % compute deriv with non-local quadrature
    Brq = DBq*DBr;
    
    % projected derivative matrix
    Dr = M\(Bq'*diag(wBq)*DBq*DBr);
    
end

% VDM for interp
BVDM = bspline_basismatrix(N+1,t,rin);


% rp = linspace(-1,1,250);
% Vp = Vandermonde1D(N,rp)/V;
% Bp = bspline_basismatrix(N+1,t,rp);
% plot(rp,Bp,'.--')
