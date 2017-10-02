clc
Globals2D

K1D = 1;
N = 3; % polynomial order


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
% return

R = getCGRestriction();

[rq sq wq] = Cubature2D(2*N);
[Vq Vrq Vsq] = bern_basis_tri(N,rq,sq); % bernstein basis
% Vq = Vandermonde2D(N,rq,sq)/V; % nodal basis
Dr = Vq\Vrq; Dr(abs(Dr)<1e-8) = 0;
Ds = Vq\Vsq; Ds(abs(Ds)<1e-8) = 0;

% mass mat
Mref = Vq'*diag(wq)*Vq;
MK = kron(spdiag(J(1,:)),Mref);
M = R*MK*R';

% % S = dudx
% DxK = kron(spdiag(rx(1,:)),Dr) + kron(spdiag(sx(1,:)),Ds);
% S = R*MK*DxK*R';
% Dx = M\S;
% Dx(abs(Dx)<1e-8) = 0;

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

cond(M)
cond(diag(1./diag(M))*M) % Jacobi
cond(diag(1./sum(M,2))*M) % lumped mass
cond(P*M)
% norm(M*P*M - M,'fro')


return

[W D] = eig(full(M));
[lam p] = sort(diag(D));
W = W(:,p);

[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = bern_basis_tri(N,rp,sp);
xp = Vp*x; yp = Vp*y;

[re se] = EquiNodes2D(N); [re se] = xytors(re,se);
Ve = bern_basis_tri(N,re,se);
xe = Ve*x; ye = Ve*y;
for i = 1:size(W,2)
    vv = reshape(R'*W(:,i),Np,K);
    clf
    color_line3(xe,ye,vv,vv,'o')
    view(3)
    pause
end
% imagesc(inv(M))

