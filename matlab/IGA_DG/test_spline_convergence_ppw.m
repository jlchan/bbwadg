clear
Globals1D;

smoothKnots = 75;'opt';

K1D = 1;

NB = 4;  Ksub = 16;
% NB = 8; Ksub = 1;

N = NB+Ksub-1; % number of sub-mesh splines

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

ndofs = K*(N+1); % total # dofs

% Initialize solver and construct grid and metric
StartUp1D;

rp = linspace(-1,1,50)';
Vp = Vandermonde1D(N,rp)/V;
xp = Vp*x;

% make splines
rB = JacobiGL(0,0,N);
xB = Vandermonde1D(N,rB)/V * x;

% [rq wq] = JacobiGQ(0,0,N);
[rgq wgq] = JacobiGQ(0,0,2*N+1);
if Ksub> 1
    [~, ~, ~, ~, ~, ~, ~, ~, VX] = bsplineVDM(NB,Ksub,rB,smoothKnots); % VDM for interp, mass, M\S
    h = @(r) repmat(diff(VX(:)'),length(r),1);
    map = @(r) reshape(h(r).*(repmat(r,1,Ksub)+1)/2 + repmat(VX(1:end-1),length(r),1),length(r)*Ksub,1);
    
    rq = map(rgq);
    wq = repmat(wgq,1,Ksub).*h(rgq)/2;
    rq = rq(:); wq = wq(:);
    %         keyboard
else
    rq = rgq(:); wq = wgq(:);
end

xq = Vandermonde1D(N,rq)/V * x;

h = 2/(K1D);
hB = h/Ksub;

% spline projection
[Bq] = bsplineVDM(NB,Ksub,rq,smoothKnots); % VDM for interp, mass, M\S
M = Bq'*spdiag(wq)*Bq;
PBq = M\(Bq'*spdiag(wq));

% polynomial projection
Vq = Vandermonde1D(NB,rq)/Vandermonde1D(NB,JacobiGL(0,0,NB));
Pq = (Vq'*spdiag(wq)*Vq)\(Vq'*spdiag(wq));

kvec = linspace(1/N,N,1000);
sk = 1;
for k = kvec
    
    %uex = @(x) cos(k*pi*x/2);
    uex = @(x) sin(k*pi*x);
    uB = PBq*uex(xq); % initialize w/lsq fit
    err = Bq*uB-uex(xq);
    err = spdiag(wq)*(err.^2).*J(1);
    L2errB(sk) = sqrt(sum(err(:)));
    sk = sk + 1;
end

% error / number of wavelengths
if smoothKnots==0
    semilogy(kvec/N,L2errB,'-','linewidth',2)
elseif strcmp(smoothKnots,'opt')==1
    semilogy(kvec/N,L2errB,'--','linewidth',2)
else
    semilogy(kvec/N,L2errB,'-.','linewidth',2)
end
hold on
set(gca,'fontsize',15);grid on
ylabel('L^2 error','fontsize',15)
axis([0 .5 10^-14 1])
return


