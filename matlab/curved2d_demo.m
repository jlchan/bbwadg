Globals2D

N = 5;
K1D = 8;
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

StartUp2D
Lx = max(VX)-min(VX);
Ly = max(VY)-min(VY);
BuildPeriodicMaps2D(Lx,Ly);

[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;

aa = .1;
x = x + aa*cos(3/2*pi*x).*sin(pi*y);
y = y + aa*sin(pi*x).*cos(3/2*pi*y);

[rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);
rxJ = rx.*J; % exactly degree N-1
sxJ = sx.*J; 
ryJ = ry.*J; 
syJ = sy.*J; 

e = ones(Nfp,1);
nrJ = [0*e; e; -e];
nsJ = [-e ; e; 0*e];

rxJf = rxJ(Fmask(:),:);
ryJf = ryJ(Fmask(:),:);
sxJf = sxJ(Fmask(:),:);
syJf = syJ(Fmask(:),:);
nxJ = diag(nrJ)*rxJf + diag(nsJ)*sxJf;
nyJ = diag(nrJ)*ryJf + diag(nsJ)*syJf;

xf = x(Fmask(:),:);
yf = y(Fmask(:),:);

% plot mesh boundaries
rp1D = linspace(-1,1,100)';
Vp1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,JacobiGL(0,0,N));
Vfp = kron(eye(Nfaces),Vp1D);
xfp = Vfp*x(Fmask(:),:);
yfp = Vfp*y(Fmask(:),:);
plot(xfp,yfp,'k.')
axis off; axis equal
set(gca,'color','none')
hold on
plot(x,y,'o')

% approximate derivative
u = .5*x.^2 + y;
dudx = x;

norm(dudx - (rxJ.*(Dr*u) + sxJ.*(Ds*u))./J,'fro')








