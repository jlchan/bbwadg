%% set up reference hex
clear
N = 5;

% interpolation nodes for geometric factors
[r1D w1D] = JacobiGL(0,0,N);
% r1D = r1D*(1-1e-10);
[r s t] = meshgrid(r1D);
r = r(:); s = s(:); t = t(:);
[wr ws wt] = meshgrid(w1D);
w = wr(:).*ws(:).*wt(:);

% % quadrature nodes for collocation
% [rq1D wq1D] = JacobiGQ(0,0,N);
% [rq sq tq] = meshgrid(rq1D);
% rq = rq(:); sq = sq(:); tq = tq(:);
% [wr ws wt] = meshgrid(wq1D);
% wq = wr(:).*ws(:).*wt(:);

% 1D matrices
V1D = Vandermonde1D(N,r1D);
D1D = GradVandermonde1D(N,r1D)/V1D;
I = eye(length(r1D));
% Vq1D = Vandermonde1D(N,rq1D)/V1D;
% M1D = Vq1D'*diag(wq1D)*Vq1D;
% Pq1D = M1D\(Vq1D'*diag(wq1D));
% Vf1D = Vandermonde1D(N,[-1;1])/Vandermonde1D(N,rq1D); % interp from quad to face points

% 3D matrices
V = kron(kron(V1D,V1D),V1D);
% Vq = kron(kron(Vq1D,Vq1D),Vq1D);
% M = Vq'*diag(wq)*Vq;
% Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));
Dr = kron(kron(I,D1D),I);
Ds = kron(kron(I,I),D1D);
Dt = kron(kron(D1D,I),I);

% % points on a quad
% [r2 s2] = meshgrid(rq1D);
% [wr2 ws2] = meshgrid(wq1D);
% r2 = r2(:); s2 = s2(:);
% e = ones(size(r2));
% wfq = wr2(:).*ws2(:);
% 
% % set face points according to r,s,t
% rf = [-e; e; r2; r2; r2; r2];
% sf = [r2; r2; -e; e; s2; s2];
% tf = [s2; s2; s2; s2; -e; e];
% wfq = repmat(wfq,6,1);

VL = Vandermonde1D(N,.5*(1+r1D)-1)/V1D;
VR = Vandermonde1D(N,.5*(1+r1D))/V1D;
Vsplit = [VL;VR];
VmL = kron(eye(N+1),VL); % interp nodes -> split mortar nodes
VmR = kron(eye(N+1),VR); % interp nodes -> split mortar nodes
Vm = [VmL;VmR];

% Vf = Vandermonde3DHex(N,rf,sf,tf)/Vandermonde3DHex(N,r,s,t); % interp nodes -> face nodes

% reference normals
nrJ = zeros((N+1)^2,6);
nsJ = zeros((N+1)^2,6);
ntJ = zeros((N+1)^2,6);
nrJ(:,1) = -1; nrJ(:,2) = 1;
nsJ(:,3) = -1; nsJ(:,4) = 1;
ntJ(:,5) = -1; ntJ(:,6) = 1;

nrJ = nrJ(:); nsJ = nsJ(:); ntJ = ntJ(:);

%% compute face indices for face points on "lines" of nodes

Fmask = zeros((N+1)^2,6);
Fmask(:,1) = find(abs(r+1)<1e-8);
Fmask(:,2) = find(abs(r-1)<1e-8);
Fmask(:,3) = find(abs(s+1)<1e-8);
Fmask(:,4) = find(abs(s-1)<1e-8);
Fmask(:,5) = find(abs(t+1)<1e-8);
Fmask(:,6) = find(abs(t-1)<1e-8);


%% make physical elements

K = 3; 
x = zeros((N+1)^3,K);
y = zeros((N+1)^3,K);
z = zeros((N+1)^3,K);

% make multi-block "mesh"
x(:,1) = .5*(1+r)-1;
y(:,1) = s;
z(:,1) = t;
x(:,2) = .5*(1+r);
% y(:,2) = s;
% z(:,2) = t;
y(:,2) = .5*(1+s)-1;
z(:,2) = t;
x(:,3) = .5*(1+r);
y(:,3) = .5*(1+s);
z(:,3) = t;

%%

% plot3(r,s,t,'o')
hold on

% e = 1;
% plot3(VmL*x(Fmask(:,2),e),VmL*y(Fmask(:,2),e),VmL*z(Fmask(:,2),e),'bo')
% text(VmL*x(Fmask(:,2),e),VmL*y(Fmask(:,2),e)-.1,VmL*z(Fmask(:,2),e),num2str((1:size(VmL,2))'))
% e = 2;
% plot3(x(Fmask(:,1),e),y(Fmask(:,1),e),z(Fmask(:,1),e),'rx')
% text(x(Fmask(:,1),e),y(Fmask(:,1),e)+.1,z(Fmask(:,1),e),num2str((1:size(VmL,2))'))

% e = 1;
% plot3(VmR*x(Fmask(:,2),e),VmR*y(Fmask(:,2),e),VmR*z(Fmask(:,2),e),'bo')
% text(VmR*x(Fmask(:,2),e),VmR*y(Fmask(:,2),e)-.1,VmR*z(Fmask(:,2),e),num2str((1:size(VmR,2))'))
% e = 3;
% plot3(x(Fmask(:,1),e),y(Fmask(:,1),e),z(Fmask(:,1),e),'rx')
% text(x(Fmask(:,1),e),y(Fmask(:,1),e)+.1,z(Fmask(:,1),e),num2str((1:size(VmR,2))'))
% return
%% curved warping?

a = 0*.1;

opt = 1;
if opt==1
    dx = cos(pi/2*x).*sin(pi*y).*sin(pi*z);
    x = x + a*dx;
    dy = sin(pi*x).*cos(pi/2*y).*sin(pi*z);
    y = y + a*dy;
    dz = sin(pi*x).*sin(pi*y).*cos(pi/2*z);
    z = z + a*dz;
elseif opt==2    
    d = (1+x).*(1-x).*(1-y).*(1+y).*(1+z).*(1-z);
    dx = d;     
    dy = d;     
    dz = 0*d;
    x = x + a*dx;
    y = y + a*dy;
    z = z + a*dz;
else    
    dx = (1-y).*(1+y).*(1+z).*(1-z);
    dy = (1+x).*(1-x).*(1+z).*(1-z);
    dz = (1+x).*(1-x).*(1-y).*(1+y);
    x = x + a*dx;
    y = y + a*dy;
    z = z + a*dz;
end

% make mesh watertight
% x(Fmask(:,1),2) = VmL*x(Fmask(:,2),1);
% y(Fmask(:,1),2) = VmL*y(Fmask(:,2),1);
% z(Fmask(:,1),2) = VmL*z(Fmask(:,2),1);
% x(Fmask(:,1),3) = VmR*x(Fmask(:,2),1);
% y(Fmask(:,1),3) = VmR*y(Fmask(:,2),1);
% z(Fmask(:,1),3) = VmR*z(Fmask(:,2),1);

%% compute geometric factors

e = ones(N+1,1); e(end) = 0; % reduce degree by 1
VNm1 = Vandermonde1D(N-1,JacobiGL(0,0,N-1));
F1D = (Vandermonde1D(N-1,JacobiGL(0,0,N))/VNm1)*(Vandermonde1D(N,JacobiGL(0,0,N-1))/V1D);
F1D = eye(N+1);

Fr = kron(kron(F1D,I),I);
Fs = kron(kron(I,F1D),I);
Ft = kron(kron(I,I),F1D);

xr = Dr*x; xs = Ds*x; xt = Dt*x;
yr = Dr*y; ys = Ds*y; yt = Dt*y;
zr = Dr*z; zs = Ds*z; zt = Dt*z;
J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);

rr = (Dr*y).*z;
rs = (Ds*y).*z;
rt = (Dt*y).*z;
% rr(Fmask(:,1),2) = VmL*rr(Fmask(:,2),1);
% rs(Fmask(:,1),2) = .5*VmL*rs(Fmask(:,2),1);
% rt(Fmask(:,1),2) = VmL*rt(Fmask(:,2),1);
% rr(Fmask(:,1),3) = VmR*rr(Fmask(:,2),1);
% rs(Fmask(:,1),3) = .5*VmR*rs(Fmask(:,2),1);
% rt(Fmask(:,1),3) = VmR*rt(Fmask(:,2),1);
% rr = Fr*rr; rs = Fs*rs; rt = Ft*rt;

rxJ = Dt*(rs) - Ds*(rt); % Q(N,N-1,N-1) + Q(N,N-1,N-1) = Q(N,N-1,N-1)
sxJ = Dr*(rt) - Dt*(rr); % Q(N-1,N,N-1)
txJ = Ds*(rr) - Dr*(rs); % Q(N-1,N-1,N)
rr = (Dr*x).*z;
rs = (Ds*x).*z;
rt = (Dt*x).*z;
% rr(Fmask(:,1),2) = VmL*rr(Fmask(:,2),1);
% rs(Fmask(:,1),2) = .5*VmL*rs(Fmask(:,2),1);
% rt(Fmask(:,1),2) = VmL*rt(Fmask(:,2),1);
% rr(Fmask(:,1),3) = VmR*rr(Fmask(:,2),1);
% rs(Fmask(:,1),3) = .5*VmR*rs(Fmask(:,2),1);
% rt(Fmask(:,1),3) = VmR*rt(Fmask(:,2),1);
% rr = Fr*rr; rs = Fs*rs; rt = Ft*rt;
ryJ = -(Dt*(rs) - Ds*(rt));
syJ = -(Dr*(rt) - Dt*(rr));
tyJ = -(Ds*(rr) - Dr*(rs));

rr = (Dr*y).*x;
rs = (Ds*y).*x;
rt = (Dt*y).*x;
% rr(Fmask(:,1),2) = VmL*rr(Fmask(:,2),1);
% rs(Fmask(:,1),2) = .5*VmL*rs(Fmask(:,2),1);
% rt(Fmask(:,1),2) = VmL*rt(Fmask(:,2),1);
% rr(Fmask(:,1),3) = VmR*rr(Fmask(:,2),1);
% rs(Fmask(:,1),3) = .5*VmR*rs(Fmask(:,2),1);
% rt(Fmask(:,1),3) = VmR*rt(Fmask(:,2),1);
% rr = Fr*rr; rs = Fs*rs; rt = Ft*rt;

rzJ = -(Dt*(rs) - Ds*(rt));
szJ = -(Dr*(rt) - Dt*(rr));
tzJ = -(Ds*(rr) - Dr*(rs));


% u = x.^2+y.^2.*z;
% dudx = 2*x;
% dudy = 2*y.*z;
% dudz = y.^2;
% dudxJ = rxJ.*(Dr*u) + sxJ.*(Ds*u) + txJ.*(Dt*u);
% dudyJ = ryJ.*(Dr*u) + syJ.*(Ds*u) + tyJ.*(Dt*u);
% dudzJ = rzJ.*(Dr*u) + szJ.*(Ds*u) + tzJ.*(Dt*u);
% norm(dudxJ - dudx.*J,'fro')
% norm(dudyJ - dudy.*J,'fro')
% norm(dudzJ - dudz.*J,'fro')

% compute normals
rxJf = rxJ(Fmask(:),:); sxJf = sxJ(Fmask(:),:); txJf = txJ(Fmask(:),:);
ryJf = ryJ(Fmask(:),:); syJf = syJ(Fmask(:),:); tyJf = tyJ(Fmask(:),:);
rzJf = rzJ(Fmask(:),:); szJf = szJ(Fmask(:),:); tzJf = tzJ(Fmask(:),:);
nxJ = nrJ.*(rxJf) + nsJ.*(sxJf) + ntJ.*(txJf);
nyJ = nrJ.*(ryJf) + nsJ.*(syJf) + ntJ.*(tyJf);
nzJ = nrJ.*(rzJf) + nsJ.*(szJf) + ntJ.*(tzJf);

fid = reshape(1:6*(N+1)^2,(N+1)^2,6);


hold on

e = 1; f = 2;
% plot3(x(:,e),y(:,e),z(:,e),'o','linewidth',2,'markersize',16)
xf1 = Vm*x(Fmask(:,f),e); 
yf1 = Vm*y(Fmask(:,f),e); 
zf1 = Vm*z(Fmask(:,f),e);
nxf1 = Vm*nxJ(fid(:,f),e); 
nyf1 = Vm*nyJ(fid(:,f),e); 
nzf1 = Vm*nzJ(fid(:,f),e);
% xf1 = x(Fmask(:,f),e); 
% yf1 = y(Fmask(:,f),e); 
% zf1 = z(Fmask(:,f),e);
% nxf1 = nxJ(fid(:,f),e); 
% nyf1 = nyJ(fid(:,f),e); 
% nzf1 = nzJ(fid(:,f),e);

plot3(xf1,yf1,zf1,'o')
quiver3(xf1,yf1,zf1,nxf1,nyf1,nzf1)

e = 2; f = 1;
% plot3(x(:,e),y(:,e),z(:,e),'x','linewidth',2,'markersize',16)
xf2 = x(Fmask(:,f),e); 
yf2 = y(Fmask(:,f),e); 
zf2 = z(Fmask(:,f),e);
nxf2 = nxJ(fid(:,f),e); 
nyf2 = nyJ(fid(:,f),e); 
nzf2 = nzJ(fid(:,f),e);
plot3(xf2,yf2,zf2,'o')
hold on
quiver3(xf2,yf2,zf2,nxf2,nyf2,nzf2)

e = 3; f = 1;
% plot3(x(:,e),y(:,e),z(:,e),'s','linewidth',2,'markersize',16)
xf3 = x(Fmask(:,f),e); 
yf3 = y(Fmask(:,f),e); 
zf3 = z(Fmask(:,f),e);
nxf3 = nxJ(fid(:,f),e); 
nyf3 = nyJ(fid(:,f),e); 
nzf3 = nzJ(fid(:,f),e);
plot3(xf3,yf3,zf3,'o')
hold on
quiver3(xf3,yf3,zf3,nxf3,nyf3,nzf3)
axis equal

norm([nxf1 + 2*[nxf2;nxf3]])
norm([nyf1 + 2*[nyf2;nyf3]])
norm([nzf1 + 2*[nzf2;nzf3]])

return

fprintf('GCL for elem  = %g\n',norm(Dr*rxJ + Ds*sxJ + Dt*txJ,'fro'))
fprintf('GCL for elem  = %g\n',norm(Dr*ryJ + Ds*syJ + Dt*tyJ,'fro'))
fprintf('GCL for elem  = %g\n',norm(Dr*rzJ + Ds*szJ + Dt*tzJ,'fro'))
