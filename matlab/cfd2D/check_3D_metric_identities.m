Globals3D
N = 3;
filename = 'Grid/cube1.msh';
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGmsh3D(filename);
VX = VX*2; VY = VY*2; VZ = VZ*2; % map -1,1
StartUp3D;

[rq sq tq wq] = tet_cubature(2*N);
Vq = Vandermonde3D(N,rq,sq,tq)/V;
M =(Vq'*diag(wq)*Vq);
Pq = M\(Vq'*diag(wq));

Drq = Vq*Dr;
Dsq = Vq*Ds;
Dtq = Vq*Dt;

Nfp = (N+1)*(N+2)/2;
[rqtri sqtri wqtri] = Cubature2D(2*N);
Vfqf = Vandermonde2D(N,rqtri,sqtri)/Vandermonde2D(N,r(Fmask(:,1)),s(Fmask(:,1)));
Vfqf = kron(eye(4),Vfqf);
rfq = Vfqf*r(Fmask(:));
sfq = Vfqf*s(Fmask(:));
tfq = Vfqf*t(Fmask(:));
wfq = repmat(wqtri,4,1);
Vfq = Vandermonde3D(N,rfq,sfq,tfq)/V;

Nfqf = length(rqtri);
e = ones(Nfqf,1); z = zeros(Nfqf,1);
nrJ = [z;z;e;-e];
nsJ = [z;-e;e;z];
ntJ = [-e;z;e;z];

[rp sp tp] = EquiNodes3D(25);
Vp = Vandermonde3D(N,rp,sp,tp)/V;

a = .1;
x = r + a*cos(pi*(r+s));
y = s - a*cos(pi*(s+t));
z = t + a*cos(pi*(r+t));

a = .05;
x = r - a*(r+s).^N;
y = s + a/2*(s+t).^N;
z = t + a/3*(r+t).^N;

% a = .1;
% x = r + a*(r+s).^floor(N/2);
% y = s + a/2*(s+t).^floor(N/2);
% z = t + a/3*(r+t).^floor(N/2);
%%
xp = Vp*x;
yp = Vp*y;
zp = Vp*z;
plot3(xp,yp,zp,'.')

% compute geofacs
xr = Dr*x; xs = Ds*x; xt = Dt*x;
yr = Dr*y; ys = Ds*y; yt = Dt*y;
zr = Dr*z; zs = Ds*z; zt = Dt*z;
J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);

% conservative curl form
rxJ = Dt*(ys.*z) - Ds*(yt.*z);
sxJ = Dr*(yt.*z) - Dt*(yr.*z);
txJ = Ds*(yr.*z) - Dr*(ys.*z);

ryJ = -(Dt*(xs.*z) - Ds*(xt.*z));
syJ = -(Dr*(xt.*z) - Dt*(xr.*z));
tyJ = -(Ds*(xr.*z) - Dr*(xs.*z));

rzJ = -(Dt*(ys.*x) - Ds*(yt.*x));
szJ = -(Dr*(yt.*x) - Dt*(yr.*x));
tzJ = -(Ds*(yr.*x) - Dr*(ys.*x));

% original cross product form
xrq = Drq*x; xsq = Dsq*x; xtq = Dtq*x;
yrq = Drq*y; ysq = Dsq*y; ytq = Dtq*y;
zrq = Drq*z; zsq = Dsq*z; ztq = Dtq*z;
rxJo =  (ysq.*ztq - zsq.*ytq); 
sxJo = -(yrq.*ztq - zrq.*ytq); 
txJo =  (yrq.*zsq - zrq.*ysq); 

ryJo = -(xsq.*ztq - zsq.*xtq); 
syJo =  (xrq.*ztq - zrq.*xtq); 
tyJo = -(xrq.*zsq - zrq.*xsq); 

rzJo = (xsq.*ytq - ysq.*xtq);
szJo = -(xrq.*ytq - yrq.*xtq);
tzJo = (xrq.*ysq - yrq.*xsq);

% div-free projection

blkdiag(L
Np = length(r);
ids = 1:Np;
Div = [Dr Ds Dt];
% Div = M\[Dr'*M Ds'*M Dt'*M];
Vv = null(Div); % find divergence free basis
wdivf = [wfq;wfq;wfq];
nhat = [nrJ;nsJ;ntJ];
Vdivfq = kron(eye(3),Vfq)*Vv; % divergence free basis evaluated on face



Vdivq = kron(eye(3),Vq)*Vv;
wdiv = [wq;wq;wq];
Mdiv = Vdivq'*diag(wdiv)*Vdivq;

PVdiv = Vv*(Mdiv \ (Vdivq'*diag(wdiv)));

rstx = PVdiv*[rxJo;sxJo;txJo]; % project onto div-free basis
rxJ = rstx(ids); sxJ = rstx(ids+Np); txJ = rstx(ids+2*Np);

rsty = PVdiv*[ryJo;syJo;tyJo]; % project onto div-free basis
ryJ = rsty(ids); syJ = rsty(ids+Np); tyJ = rsty(ids+2*Np);

rstz = PVdiv*[rzJo;szJo;tzJo]; % project onto div-free basis
rzJ = rstz(ids); szJ = rstz(ids+Np); tzJ = rstz(ids+2*Np);

sqrt(sum(wq.*(Vq*rxJ-rxJo).^2))
sqrt(sum(wq.*(Vq*ryJ-ryJo).^2))
sqrt(sum(wq.*(Vq*rzJ-rzJo).^2))
% norm(rxJ-Pq*rxJo)+norm(sxJ-Pq*sxJo)+norm(txJ-Pq*txJo)+norm(ryJ-Pq*ryJo)+norm(syJ-Pq*syJo)+norm(tyJ-Pq*tyJo)+norm(rzJ-Pq*rzJo)+norm(szJ-Pq*szJo)+norm(tzJ-Pq*tzJo)

% return
% check metric identities
norm(Dr*rxJ + Ds*sxJ + Dt*txJ)
norm(Dr*ryJ + Ds*syJ + Dt*tyJ)
norm(Dr*rzJ + Ds*szJ + Dt*tzJ)



