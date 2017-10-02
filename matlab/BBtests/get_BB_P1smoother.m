function W = get_BB_P1smoother(Nin)

% Nin = 5;
N = Nin;

K1D = N; % macro grid for P1 FEM
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

xq = Vq*x;

V_sub = bern_basis_1D(N,xq(:)); % interpolate to micro grid
wJ = diag(wq)*(Vq*J);

% P1-FEM based kantorovich
V1 = Vandermonde1D(1,rq)/Vandermonde1D(1,[-1;1]);
V1 = kron(eye(N),V1);
R{1} = 1;
for i = 2:N
    R{i} = [1;1];
end
R{N+1} = 1;
R = blkdiag(R{:});

W = (V1*R)'*diag(wJ(:))*V_sub * (N+1)/2;
W = W';
% sum(W,1)
% sum(W,2)
% keyboard



