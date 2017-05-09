Globals2D

N = 2;
K1D = 1;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

% Initialize solver and construct grid and metric
StartUp2D;

a = 3/8;
x = x + a*cos(.5*3*pi*y).*cos(pi/2*x);
y = y + a*sin(.5*3*pi*x).*cos(pi/2*y);

[rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);
[nx, ny, sJ] = Normals2D();

Nq = 2*N;
[rq sq wq] = Cubature2D(Nq);
Vq = Vandermonde2D(N,rq,sq)/V;
[Vrq Vsq] = GradVandermonde2D(N,rq,sq);
Vrq = Vrq/V; Vsq = Vsq/V;

[rq1D wq1D] = JacobiGQ(0,0,Nq);
wfq = [wq1D;wq1D;wq1D];
rfq = [rq1D; -rq1D; -ones(size(rq1D))];
sfq = [-ones(size(rq1D));rq1D; -rq1D];
% plot(rfq,sfq,'o');return

Vfq = Vandermonde2D(N,rfq,sfq)/V;

Vfqf = Vandermonde1D(N,rq1D)/Vandermonde1D(N,JacobiGL(0,0,N));
Vfqf = blkdiag(Vfqf,Vfqf,Vfqf);

norm(Dr*(rx.*J) + Ds*(sx.*J),'fro')
norm(Dr*(ry.*J) + Ds*(sy.*J),'fro')

rxJ = rx.*J; sxJ = sx.*J;
ryJ = ry.*J; syJ = sy.*J;

rxJ = rxJ(Fmask(:),:); sxJ = sxJ(Fmask(:),:);
ryJ = ryJ(Fmask(:),:); syJ = syJ(Fmask(:),:);

nxJ = nx.*sJ;
nyJ = ny.*sJ;

return

rxJ = Vq*(rx.*J); sxJ = Vq*(sx.*J);
ryJ = Vq*(ry.*J); syJ = Vq*(sy.*J);
J = Vq*J;
nxJ = Vfqf*(nx.*sJ); 
nyJ = Vfqf*(ny.*sJ);

e = 1;
u = x(:,e).^N + y(:,e).^N;
dudxJ = rxJ(:,e).*(Vrq*u) + sxJ(:,e).*(Vsq*u);

r1 = Vq'*diag(wq)*dudxJ;
r2 = Vrq'*diag(rxJ(:,e).*wq)*Vq*u + Vsq'*diag(sxJ(:,e).*wq)*Vq*u - Vfq'*diag(wfq.*nxJ(:,e))*Vfq*u;
norm(r1 + r2)
