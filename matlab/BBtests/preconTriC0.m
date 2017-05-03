% clc
clear
Globals2D

K1D = 4;

% filename = 'Grid/Other/block2.neu'; 
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename); % unstructured grid
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D,K1D); % structured grid

%% make P1 precon
if 1
    N = 1;
    StartUp2D;
    R1 = getCGRestriction();
    
    V1 = V;
    Mref1 = inv(V1*V1');
    MK1 = kron(spdiag(J(1,:)),Mref1);
    M1 = R1*MK1*R1'; % make low order mass matrix
    invM1 = sum(1./diag(sum(M1,2))); 
end

%% 

N = 9; % polynomial order

StartUp2D;
PlotMesh2D

[re se] = EquiNodes2D(N); [re se] = xytors(re,se);
Ve = Vandermonde2D(N,re,se)/V;
xe = Ve*x; ye = Ve*y;

[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

R = getCGRestriction();
[gids,~] = find(R);
gids = reshape(gids,Np,K);

[rq sq wq] = Cubature2D(2*N);
% Vq = bern_basis_tri(N,rq,sq); % bernstein basis
Vq = Vandermonde2D(N,rq,sq)/V; % nodal basis
Veq = Vandermonde2D(N,rq,sq)/Ve; % nodal basis

% mass mat
Mref = Vq'*diag(wq)*Vq;
MK = kron(spdiag(J(1,:)),Mref);
M = R*MK*R';

[r1 s1] = Nodes2D(1); [r1 s1] = xytors(r1,s1);
V1 = Vandermonde2D(1,r1,s1);
VDM1 = Vandermonde2D(1,r,s)/V1;
E1N = kron(speye(K),VDM1); % extend from P1 to PN
% P1 = R*kron(speye(K),V1);

% P1 lumped correction - 
PK = kron(spdiag(1./J(1,:)),inv(Mref));
P = diag(1./sum(R,2))*R* ( E1N*R1'*(diag(1./diag(M1))*(R1*E1N'*MK)) + PK ) * R'; % project onto P1
P = (M1\(R1*E1N'*MK*R'))'*diag(1./sum(M1,2))*(M1\(R1*E1N'*MK*R')) + diag(1./sqrt(sum(R,2)))*R*PK*R'* diag(1./sqrt(sum(R,2))); % project onto P1

cond(diag(1./diag(M))*M) % Jacobi
cond(diag(1./sum(M,2))*M) % mass lump
cond(P*M)
% cond(M)

return



invM = inv(M);
for i = 1:size(invM,2)
    vv = Vp*reshape(R'*invM(:,i),Np,K);
    clf
color_line3(xp,yp,vv,vv,'.')
view(3)
axis normal
pause
end