function checkNodalPyr

clear
N = 4;

[rq sq tq w] = pyr_cubature(N);
[Vq] = pyr_nodal_basis(N,rq,sq,tq);

[V] = pyr_basis(N,rq,sq,tq);
norm(V - Vq*(Vq\V),'fro')

% mapped

aa = .2;  ids = [1 3 4 2 5];
[VX VY VZ] = pyr_nodes(1); VX = VX(ids); VY = VY(ids); VZ = VZ(ids);
VX = VX + aa*randn(size(VX)); VY = VY + aa*randn(size(VX)); VZ = VZ + aa*randn(size(VX));
[xq,yq,zq,rx,sx,tx,ry,sy,ty,rz,sz,tz,J] = pyr_geom_factors(VX,VY,VZ,rq,sq,tq);

wJ = w.*J;
M = Vq'*diag(wJ)*Vq;
M(abs(M)<1e-8)=0;
spy(M)

% orthogonal basis on mapped pyramidal elements. semi nodal

function [V] = pyr_nodal_basis(N,r,s,t)

% convert to abc
a = 2*(r+1)./(1-t)-1;
b = 2*(s+1)./(1-t)-1;
c = t;

ids = abs(1-t)<1e-12;
a(ids) = -1; b(ids) = -1;

% change of vars from a to b
dadr = 2./(1-t);
dbds = 2./(1-t);
dadt = (2*r + 2)./(1 - t).^2;
dbdt = (2*s + 2)./(1 - t).^2;
dcdt = 1;

Np = (N+1)*(N+2)*(2*N+3)/6;
V = zeros(length(r),Np);
Vr = V; Vs = V; Vt = V;

ind = 1;
for k = 0:N
    [bq aq wab] = QCubature2D(k);
    VDM = QVandermonde2D(k,aq,bq);
    [Vab Va Vb] = QVandermonde2D(k,a,b);
    
    Vab = Vab/VDM;
    
    CNk = (N+2)/(2^(2*k+2)*(2*k+3));
%     pNk = JacobiP(c,2*k+3,0,N-k);
%     pc = ((1-c)/2).^k.*pNk;
    c1D = JacobiGQ(1,0,N); invVc = inv(Vandermonde1D(N,c1D));
    pc = Vandermonde1D(N,c)*invVc(:,k+1);
    
    for ij = 1:(k+1)^2
        scale = 1/sqrt(CNk*wab(ij)); % normalization
        V(:,ind) = scale*Vab(:,ij).*pc;
        ind = ind+1;
    end
end

