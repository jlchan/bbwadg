function W = get_BB_smoother(Nin)
% Nin = 3;
N = Nin;

K1D = N+1; % macro grid
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

r = JacobiGL(0,0,N);
va = EToV(:,1)'; vb = EToV(:,2)';
x = ones(N+1,1)*VX(va) + 0.5*(r+1)*(VX(vb)-VX(va));

V  = Vandermonde1D(N, r); 
Dr = GradVandermonde1D(N,r)/V;

[rx,J] = GeometricFactors1D(x,Dr); % calculate geometric factors

[rq wq] = JacobiGQ(0,0,N+1);
Vq = Vandermonde1D(N,rq)/V;

re = linspace(-1,1,N+1)';
rp = linspace(-1,1,200)';
VB = bern_basis_1D(N,r);
VBe = bern_basis_1D(N,re);

V_sub = bern_basis_1D(N,x(:)); % interpolate to micro grid
wJ = diag(wq)*(Vq*J);

utmp = zeros(N+1,1);
W = zeros(N+1);
for i = 1:N+1
    utmp(i) = 1;
    uB = (N+1)/2 * sum(wJ.*(Vq*reshape(V_sub*utmp,N+1,K)),1);
    W(:,i) = uB(:);
    utmp(i) = 0;
end

% preserve linears exactly
if 0
    [V D] = eig(W);
    V = real(V);
    lam = real(diag(D));
    [~, p] = sort(lam);
    lam(p(N)) = 1; % try to preserve linears
    W = V*diag(lam)/V;
end

