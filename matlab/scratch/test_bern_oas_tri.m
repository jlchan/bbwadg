Globals2D

% Polynomial order used for approximation 
N = 3;

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/squarereg.neu');
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/squareireg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/block2.neu');

% Initialize solver and construct grid and metric
StartUp2D;

% map plotting nodes
[rp sp] = EquiNodes2D(75); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

% cubature for local mass
[rq sq w] = Cubature2D(2*N);

% switch b/w Bern and W&B nodal 
if 1
    [V Vr Vs V1 V2 V3 id] = bern_basis_tri(N,r,s);
    Vq = bern_basis_tri(N,rq,sq);
    Vp = bern_basis_tri(N,rp,sp);
    Dr = V\Vr; Dr(abs(Dr)<1e-8) = 0;
    Ds = V\Vs; Ds(abs(Ds)<1e-8) = 0;
else   
    Vq = Vandermonde2D(N,rq,sq)/V;    
    Vp = Vandermonde2D(N,rp,sp)/V;    
end

% greville abscissae
[rg sg] = EquiNodes2D(N); [rg sg] = xytors(rg,sg);

% for i = 1:Np
%     vv = Vp(:,i);
%     clf;   hold on
%     color_line3(rp,sp,vv,vv,'.');
%     plot(rg,sg,'go','markersize',32)
%     text(rg+.05,sg,num2str((1:Np)'))
%     pause
% end

[R vmapBT] = getCGRestriction(); R = R';

% get element-constant geofacs
JK = J(1,:);
rxK = rx(1,:); sxK = sx(1,:);
ryK = ry(1,:); syK = sy(1,:);

% get block unassembled operators
MK = Vq'*diag(w)*Vq;
M = kron(spdiag(JK),MK); 
Dx = kron(spdiag(rxK),Dr) + kron(spdiag(sxK),Ds);
Dy = kron(spdiag(ryK),Dr) + kron(spdiag(syK),Ds);
Ks = Dx'*M*Dx + Dy'*M*Dy;

% forcing = 1
b = sum(M,2);

A = R'*Ks*R;
b = R'*b;

keyboard
% b = b*0;
% b(25) = 1;



b(vmapBT) = 0;
A(vmapBT,:) = 0; A(:,vmapBT) = 0;
A(vmapBT,vmapBT) = speye(length(vmapBT));

u = reshape(R*(A\b),Np,K);

vv = Vp*u;
figure
color_line3(xp,yp,vv,vv,'.')