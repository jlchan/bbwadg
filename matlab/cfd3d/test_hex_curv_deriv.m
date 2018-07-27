% clear
N = 2;
[r1D w1D] = JacobiGL(0,0,N);
[r s t] = meshgrid(r1D);
r = r(:); s = s(:); t = t(:);

NODETOL = 1e-8;

[r2 s2] = meshgrid(r1D);
r2 = r2(:);
s2 = s2(:);
e = ones(size(r2));

rf = [-e; e; r2; r2; r2; r2];
sf = [r2; r2; -e; e; s2; s2];
tf = [s2; s2; s2; s2; -e; e];

Vf = Vandermonde3DHex(N,rf,sf,tf)/Vandermonde3DHex(N,r,s,t);
NODETOL = 1e-8;

x = r;
y = s;
z = t;

a = .1;
d = exp((x+y+z)/8);
x = x + a*d;
y = y + a*d;
z = z + a*d;

V1D = Vandermonde1D(N,r1D);
D1D = GradVandermonde1D(N,r1D)/V1D;
I = eye(length(r1D));
Dr = kron(kron(I,D1D),I);
Ds = kron(kron(I,I),D1D);
Dt = kron(kron(D1D,I),I);

xr = Dr*x; xs = Ds*x; xt = Dt*x;
yr = Dr*y; ys = Ds*y; yt = Dt*y;
zr = Dr*z; zs = Ds*z; zt = Dt*z;
J = xr.*(ys.*zt-zs.*yt) - ...
    yr.*(xs.*zt-zs.*xt) + ...
    zr.*(xs.*yt-ys.*xt);

rxJ = Dt*((Ds*y).*z) - Ds*((Dt*y).*z); % this is the problematic one
sxJ = Dr*((Dt*y).*z) - Dt*((Dr*y).*z);
txJ = Ds*((Dr*y).*z) - Dr*((Ds*y).*z);

ryJ = -(Dt*((Ds*x).*z) - Ds*((Dt*x).*z));
syJ = -(Dr*((Dt*x).*z) - Dt*((Dr*x).*z));
tyJ = -(Ds*((Dr*x).*z) - Dr*((Ds*x).*z));

rzJ = -(Dt*((Ds*y).*x) - Ds*((Dt*y).*x));
szJ = -(Dr*((Dt*y).*x) - Dt*((Dr*y).*x));
tzJ = -(Ds*((Dr*y).*x) - Dr*((Ds*y).*x));

fprintf('GCL for elem  = %g\n',norm(Dr*rxJ + Ds*sxJ + Dt*txJ,'fro'))

u = x;
duex = 1.0 + 0*x;

d1 = rxJ.*(Dr*u) + sxJ.*(Ds*u) + txJ.*(Dt*u);
d2 = Dr*(rxJ.*u) + Ds*(sxJ.*u) + Dt*(txJ.*u);
dudx = .5*(d1+d2)./J;

sum(sum(abs(dudx-duex)))


