Globals2D

N = 3;
K1D = 1;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
StartUp2D;

% [r s] = Nodes2D(N); [r s] = xytors(r,s);
% V = Vandermonde2D(N,r,s);
% [Vr Vs] = GradVandermonde2D(N,r,s);
% Dr = Vr/V; Ds = Vs/V;

[rq sq wq] = Cubature2D(2*N+1);
Vq = Vandermonde2D(N,rq,sq)/V;

M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq));
W = diag(wq);

[rq1D wq1D] = JacobiGQ(0,0,N);
rfq = [rq1D; -rq1D; -ones(size(rq1D))];
sfq = [-ones(size(rq1D)); rq1D; -rq1D];
wfq = [wq1D; wq1D; wq1D];
Vq1D = Vandermonde1D(N,rq1D)/Vandermonde1D(N,JacobiGL(0,0,N));
% plot(rfq,sfq,'o')

Vfq = Vandermonde2D(N,rfq,sfq)/V;
Vfqf = kron(eye(3),Vq1D);
Mf = Vfq'*diag(wfq)*Vfq;
L = M\(Vfq'*diag(wfq));
% L = M\Mf;

%nr = [-zeros(size(rq1D)); ones(size(rq1D))/sqrt(2); -ones(size(rq1D)); ];
%ns = [-ones(size(rq1D)); ones(size(rq1D))/sqrt(2); -zeros(size(rq1D)); ];
nrJ = [-zeros(size(rq1D)); ones(size(rq1D)); -ones(size(rq1D)); ];
nsJ = [-ones(size(rq1D)); ones(size(rq1D)); -zeros(size(rq1D)); ];
% quiver(rfq,sfq,nr,ns)
% hold on
% plot(rfq,sfq,'o')

norm(W*Vq*Dr*Pq - (W*Vq*L*diag(nrJ)*Vfq*Pq - Pq'*Dr'*Vq'*W ),'fro')
norm(W*Vq*Ds*Pq - (W*Vq*L*diag(nsJ)*Vfq*Pq - Pq'*Ds'*Vq'*W ),'fro')

xfq = Vfq*x; yfq = Vfq*y;
nxq = Vfqf*nx; nyq = Vfqf*ny;
% plot(xfq,yfq,'o')
% hold on
% quiver(xfq,yfq,nxq,nyq)

rxJ = rx.*J; sxJ = sx.*J;
ryJ = ry.*J; syJ = sy.*J;

rxJf = Vfq*rxJ; sxJf = Vfq*sxJ;
ryJf = Vfq*ryJ; syJf = Vfq*syJ;

nxJ = rxJf.*repmat(nrJ,1,K) + sxJf.*repmat(nsJ,1,K);
nyJ = ryJf.*repmat(nrJ,1,K) + syJf.*repmat(nsJ,1,K);
sJq = Vfqf*sJ;

% clf
% plot(xfq,yfq,'o')
% hold on
% quiver(xfq,yfq,nxJ,nyJ)
 
norm(Vfqf*(nx.*sJ) - nxJ,'fro')
norm(Vfqf*(ny.*sJ) - nyJ,'fro')


