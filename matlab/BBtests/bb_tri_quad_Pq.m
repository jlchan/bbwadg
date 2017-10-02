clear

N = 2;
[r s w a b] = tri_tp_cubature(N);

V = bern_basis_tri(N,r,s);
M = V'*diag(w)*V;

[a1D w1D] = JacobiGQ(0,0,N);

VBquad = zeros(length(a),(N+1)^2);
Va = bern_basis_1D(N,a1D);

[Eth Etw Etq] = get_Eth(N);

sk = 1;
for j = 0:N
    for i = 0:N
        Va = bern_basis_1D(N,a);
        Vb = bern_basis_1D(N,b);
        VBquad(:,sk) = Va(:,i+1).*Vb(:,j+1);
        sk = sk + 1;
    end
end

norm(V-VBquad*Etq,'fro')

% steps to use Pq below:
% - Etq: elevate tri to quad
% - TP quadrature + project onto TP quad BB basis.
% - integrate result against tri basis using mass matrix: Etq'*Mquad*u = (BNtri_i,u)
% pros: well conditioned. 
% cons: can't figure out fast way to apply

% TO CHECK: top block row: 1D degree reduction = lower block rows.

% [wa wb] = meshgrid(w1D); wquad = wa(:).*wb(:);
wquad = w; % project wrt triangle Jacobi weight
Mquad = VBquad'*diag(wquad)*VBquad; 

Pq = M\(Etq'*Mquad);
Pq(abs(Pq)<1e-8) = 0;

for i = 0:N
    E{i+1} = bern_basis_1D(N,a1D)\bern_basis_1D(N-i,a1D);
    if i < N
        Ei{i+1} = bern_basis_1D(N-i+1,a1D)\bern_basis_1D(N-i,a1D);
    end
end
U = Pq - pinv(Etq);
U(abs(U)<1e-8) = 0;

T = bern_basis_tri(N,r,s)\Vandermonde2D(N,r,s);
T1D = bern_basis_1D(N,a1D)\Vandermonde1D(N,a1D);
% cond(Pq)


ids = 1:N+1;
for i = 1:N
    off = i*(N+1);
    Ui{i} = U(ids,ids+off);
end
% E{2}'*U(ids,ids+off)
iEtq = pinv(Etq);