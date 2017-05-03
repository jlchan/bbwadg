clc
clear
Globals2D

K1D = 1;
N = 2; % polynomial order

% filename = 'Grid/Other/block2.neu'; 
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename); % unstructured grid
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D,K1D); % structured grid

StartUp2D;

[re se] = EquiNodes2D(N); [re se] = xytors(re,se);
Ve = Vandermonde2D(N,re,se)/V;
x = Ve*x; y = Ve*y;

[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;


PlotMesh2D
hold on
plot(x,y,'o')
off = zeros(Np,K);
off(:,1) = .05;
off(:,2) = -.05;
text(x(:)+.05,y(:)+off(:),num2str((1:Np*K)'))

ids(1) = 1;
for i = 1:N
    ids(i+1) = ids(i) + (N+1-(i-1));
end
ids1 = ids(1:end-1)+1;
ids2 = ids(1:end-1)+1 + Np;

vol_ids = 1:Np;
%I1 = [vol_ids ids2];
%I2 = [vol_ids+Np ids1];

I1 = [vol_ids];
I2 = [vol_ids+Np];
%plot(x(I1),y(I1),'rs','markersize',10)
plot(x(I2),y(I2),'rs','markersize',12)

R = getCGRestriction();

[rq sq wq] = Cubature2D(2*N);
Vq = bern_basis_tri(N,rq,sq); % bernstein basis
% Vq = Vandermonde2D(N,rq,sq)/V; % nodal basis

% mass mat
Mref = Vq'*diag(wq)*Vq;
MK = kron(spdiag(J(1,:)),Mref);
M = R*MK*R';

VDMe = kron(eye(K),Ve)*R';

vv = VDMe*rand(size(R',2),1);
vv = Vp*reshape(vv,Np,K);
color_line3(xp,yp,vv,vv,'.')

return

% precon
PK = kron(spdiag(1./J(1,:)),inv(Mref));
P = R*PK*R';
P = diag(1./sum(R,2))*R*PK*R'; % scaled precon

P = zeros(size(M));
[gids,~] = find(R);
gI1 = gids(I1);
gI2 = gids(I2);
P(gI1,gI1) = P(gI1,gI1) + inv(M(gI1,gI1));
P(gI2,gI2) = P(gI2,gI2) + inv(M(gI2,gI2));

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

% cond(diag(1./diag(M))*M) % Jacobi
% cond(diag(1./sum(M,2))*M) % lumped mass
cond(P*M)
% norm(M*P*M - M,'fro')

% imagesc(inv(M))

