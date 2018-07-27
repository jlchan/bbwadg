function [LIFT] = LiftQuad2D_DGSEM()

% function [LIFT] = LiftQuad2D()
% Purpose  : Compute surface to volume lift term for DG formulation

Globals2D;
Emat = zeros(Np, Nfaces*Nfp);

% face 1
faceR = r(Fmask(:,1));
V1D = Vandermonde1D(N, faceR); 
massEdge1 = inv(V1D*V1D');
Emat(Fmask(:,1),1:Nfp) = diag(sum(massEdge1,2));

% face 2
faceR = s(Fmask(:,2));
V1D = Vandermonde1D(N, faceR);
massEdge2 = inv(V1D*V1D');
Emat(Fmask(:,2),Nfp+1:2*Nfp) = diag(sum(massEdge2,2));

% face 3
faceR = r(Fmask(:,3));
V1D = Vandermonde1D(N, faceR); 
massEdge3 = inv(V1D*V1D');
Emat(Fmask(:,3),2*Nfp+1:3*Nfp) = diag(sum(massEdge3,2));

% face 4
faceS = s(Fmask(:,4));
V1D = Vandermonde1D(N, faceS); 
massEdge4 = inv(V1D*V1D');
Emat(Fmask(:,4),3*Nfp+1:4*Nfp) = diag(sum(massEdge4,2));

% inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
M = inv(V*V');
LIFT = diag(1./sum(M,2))*Emat;

ids = find(abs(LIFT)<1e-13); LIFT(ids) = 0;
LIFT = sparse(LIFT);

return