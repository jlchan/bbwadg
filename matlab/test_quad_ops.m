clear
Globals2D

N = 2;
[Nv, VX, VY, K, EToV] = QuadMesh2D(2,2);

StartUp2D

% x = x + .1*cos(pi/2*x).*cos(3/2*pi*y);
% y = y + .1*cos(3/2*pi*x).*cos(pi/2*y);
% [rx,sx,ry,sy,J] = GeometricFactorsQuad2D(x,y,Dr,Ds);
% [nx, ny, sJ] = NormalsQuad2D();

[rq1D wq1D] = JacobiGQ(0,0,N);
e = ones(size(rq1D));
rf = [rq1D;e;-rq1D;-e];
sf = [-e;rq1D;e;-rq1D];
wf = [wq1D;wq1D;wq1D;wq1D];
nrJ = [0*e;e;0*e;-e];
nsJ = [-e;0*e;e;0*e];
Vf = Vandermonde2D(N,rf,sf)/V;
xf = Vf*x;
yf = Vf*y;

D1D = GradVandermonde1D(N,rq1D)/Vandermonde1D(N,rq1D);

V1D = Vandermonde1D(N,rq1D)/Vandermonde1D(N,JacobiGL(0,0,N));
Vface = kron(eye(Nfaces),V1D);
Vf1D = Vandermonde1D(N,[-1,1])/Vandermonde1D(N,rq1D);

[rq sq] = meshgrid(rq1D); rq = rq(:); sq = sq(:);
Vq = Vandermonde2D(N,rq,sq)/V;
% kron(D1D,eye(N+1))
% kron(eye(N+1), D1D)


Vf = Vandermonde2D(N,rf,sf)/Vandermonde2D(N,rq,sq);
[wr ws]= meshgrid(wq1D);
wq = wr(:).*ws(:);
Lf = diag(1./wq)*Vf'*diag(wf);
Lf(abs(Lf)<1e-8) = 0;
Lf1D = diag(1./wq1D)*Vf1D';