clear
N = 4;

[r s] = Nodes2D(N); [r s] = xytors(r,s);
%[rq sq wq] = Cubature2D(2*N+1);
[rq sq wq a b] = tri_tp_cubature(2*N+1);

V = Vandermonde2D(N,r,s);
Vq = Vandermonde2D(N,rq,sq);
norm(Vq'*diag(wq)*Vq - eye(size(Vq,2)),'fro')

VB = bern_basis_tri(N,r,s);
VBq = bern_basis_tri(N,rq,sq);

M = VBq'*diag(wq)*VBq;

% modal BB coefficients - orthog wrt l2
T = VB\V;
T(abs(T)<1e-8) = 0;

A = T'*T;
A(abs(A)<1e-8) = 0;

[re se] = EquiNodes2D(N); [re se] = xytors(re,se);
t = T(:,6);
plot3(re,se,t,'.','markersize',32)

[Eth Etw Etq] = get_Eth(N);

Tq = Etq*T; % modal tri to bb quad
[ae be] = meshgrid(linspace(-1,1,N+1));
ae = ae(:); be = be(:)

t = Tq(:,4);
plot3(ae,be,t,'.','markersize',32)