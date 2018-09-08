% clear
N = 4;
[r1D w1D] = JacobiGL(0,0,N);
[r s t] = meshgrid(r1D);
r = r(:); s = s(:); t = t(:);
[wr ws wt] = meshgrid(w1D);
w = wr(:).*ws(:).*wt(:);

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

a = .25;
d = exp((x+y+z)/8);
x = x + a*d;
y = y + a*d;
z = z + a*d;

V1D = Vandermonde1D(N,r1D);
V = kron(kron(V1D,V1D),V1D);
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

Fr = (Dr*y).*z;
Fs = (Ds*y).*z;
Ft = (Dt*y).*z;
e = ones(N+1,1); e(end) = 0; % reduce degree by 1
VNm1 = Vandermonde1D(N-1,JacobiGL(0,0,N-1));
F1D = (Vandermonde1D(N-1,JacobiGL(0,0,N))/VNm1)*(Vandermonde1D(N,JacobiGL(0,0,N-1))/V1D); 

Fr = kron(kron(I,F1D),I)*Fr;
Fs = kron(kron(I,I),F1D)*Fs;
Ft = kron(kron(F1D,I),I)*Ft;
rxJ = Dt*(Fs) - Ds*(Ft); 
sxJ = Dr*(Ft) - Dt*(Fr);
txJ = Ds*(Fr) - Dr*(Fs);

ryJ = -(Dt*((Ds*x).*z) - Ds*((Dt*x).*z));
syJ = -(Dr*((Dt*x).*z) - Dt*((Dr*x).*z));
tyJ = -(Ds*((Dr*x).*z) - Dr*((Ds*x).*z));

rzJ = -(Dt*((Ds*y).*x) - Ds*((Dt*y).*x));
szJ = -(Dr*((Dt*y).*x) - Dt*((Dr*y).*x));
tzJ = -(Ds*((Dr*y).*x) - Dr*((Ds*y).*x));

% u = x;
% duex = 1.0 + 0*x;
% d1 = rxJ.*(Dr*u) + sxJ.*(Ds*u) + txJ.*(Dt*u);
% d2 = Dr*(rxJ.*u) + Ds*(sxJ.*u) + Dt*(txJ.*u);
% dudx = .5*(d1+d2)./J;
% sum(sum(abs(dudx-duex)))
% fprintf('GCL for elem  = %g\n',norm(Dr*rxJ + Ds*sxJ + Dt*txJ,'fro'))

u = ((1+r).*(.5+s).*(.25+t)).^N;
v = (2+r).^N.*((1+s).*(.5+t)).^(N-1);
v = rxJ;

%i1 = sum(w.*(rxJ.*(Dr*u)));
i1 = sum(w.*(v.*(Dr*u)));

[rq1D wq1D] = JacobiGQ(0,0,N+1);
[rq sq tq] = meshgrid(rq1D);
rq = rq(:); sq = sq(:); tq = tq(:);
[wr ws wt] = meshgrid(wq1D);
wq = wr(:).*ws(:).*wt(:); 
Vq1D = Vandermonde1D(N,rq1D)/V1D;
Vq = kron(kron(Vq1D,Vq1D),Vq1D);
i2 = sum(wq.*((Vq*v).*(Vq*Dr*u)));

i1-i2


