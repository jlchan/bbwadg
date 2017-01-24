N = 6;
[r s] = Nodes2D(N); [r s] = xytors(r,s);
Np = length(r);
[rq sq wq] = Cubature2D(2*N);

Vq = bern_basis_tri(N,rq,sq);
% Vq = Vandermonde2D(N,rq,sq)/Vandermonde2D(N,r,s);
M = Vq'*diag(wq)*Vq;

plot(r,s,'o'); text(r,s,num2str((1:length(r))'))

vids = [1; N+1; Np];
fids = unique([
    find(abs(s+1)<1e-8);
    find(abs(r+1)<1e-8); 
    find(abs(r+s)<1e-8)]);

fids = setdiff(fids,vids);
iids = setdiff(1:Np,[vids;fids]);

gids = [fids;vids];
gids = sort(gids);

A = M(iids,iids);
B = M(iids,gids);
C = M(gids,gids);

S = C - B'*(A\B);

r1D = JacobiGL(0,0,N);
[r1Dq w1Dq] = JacobiGL(0,0,N);
L = bern_basis_1D(N,r1D)\Vandermonde1D(N,r1D);
V1Dq = bern_basis_1D(N,r1Dq);
M = V1Dq'*diag(w1Dq)*V1Dq;