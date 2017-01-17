clc
Globals2D

K1D = 1;
N = 4; % polynomial order


% filename = 'Grid/Other/block2.neu'; 
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename); % unstructured grid
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D,K1D); % structured grid

StartUp2D;

PlotMesh2D
hold on
plot(x,y,'o')
off = zeros(Np,K);
off(:,1) = .05;
off(:,2) = -.05;
text(x(:)+.05,y(:)+off(:),num2str((1:Np*K)'))
return

R = getCGRestriction();

[rq sq wq] = Cubature2D(2*N);
% Vq = bern_basis_tri(N,rq,sq); % bernstein basis
Vq = Vandermonde2D(N,rq,sq)/V; % nodal basis

% mass mat
Mref = Vq'*diag(wq)*Vq;
MK = kron(spdiag(J(1,:)),Mref);
M = R*MK*R';

% precon
PK = kron(spdiag(1./J(1,:)),inv(Mref));
P = R*PK*R';
P = diag(1./sum(R,2))*R*PK*R'; % scaled precon

if 1 % check mass matrix reproduces global polynoms
    f = x.^N+y.^N;
    u = M\(R*MK*f(:));
    u = reshape(R'*u,Np,K);
    err = norm(u-f,'fro');
    if err > 1e-10
        err
        keyboard
    end
end

cond(diag(1./diag(M))*M) % Jacobi
cond(diag(1./sum(M,2))*M) % lumped mass
cond(P*M)
% norm(M*P*M - M,'fro')

% imagesc(inv(M))

