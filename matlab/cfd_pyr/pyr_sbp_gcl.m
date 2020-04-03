clear
N = 5;

[r s t] = pyr_nodes(N);
t = t*(1-1e-8);
[rq sq tq wq] = pyr_cubature(2*N);
% [rp sp tp] = pyr_equi_nodes(20);
[ap bp cp] = meshgrid(linspace(-1,1-2/N,20));
% ap = ap(:); bp = bp(:); cp = cp(:);
[rp sp tp] = pyr_abctorst(ap,bp,cp)

[rt st wt] = Cubature2D(2*N);
et = ones(size(rt));
zt = 0*et;
[rq1D wq1D] = JacobiGQ(0,0,N);
[rquad squad] = meshgrid(rq1D);
rquad = rquad(:); 
squad = squad(:);
equad = ones(size(rquad));
zquad = 0*equad;
[wr ws] = meshgrid(wq1D); 
wquad = wr(:).*ws(:);

rf = [ rt; -rt;  rt; -et;  rquad];
sf = [-et;  st; -st;  rt;  squad];
tf = [ st;  rt;  st;  st; -equad];
wf = [ wt;  wt;  wt;  wt;  wquad];

nrJq = [ zt; et; zt; -et;  zquad];
nsJq = [-et; zt; et;  zt;  zquad];
ntJq = [ zt; et; et;  zt; -equad];

tol = 1e-7;
fids1 = find(abs(s+1)<tol);
fids2 = find(abs(r+t)<tol);
fids3 = find(abs(s+t)<tol);
fids4 = find(abs(r+1)<tol);
fids5 = find(abs(t+1)<tol);
Fmask = [fids1; fids2; fids3; fids4; fids5];
et = ones(size(fids1));
zt = 0*et;
nrJ = [ zt; et; zt; -et;  zquad];
nsJ = [-et; zt; et;  zt;  zquad];
ntJ = [ zt; et; et;  zt; -equad];

% [Vq Vrq Vsq Vtq] = bern_basis_pyr(N,rq,sq,tq);
[V Vr Vs Vt] = pyr_basis(N,r,s,t);
[Vq Vrq Vsq Vtq] = pyr_basis(N,rq,sq,tq);
Vq = Vq/V;
Vrq = Vrq/V;
Vsq = Vsq/V;
Vtq = Vtq/V;

[Vp Vrp Vsp Vtp] = pyr_basis(N,rp,sp,tp);
Vp = Vp/V;
Drp = Vrp/V;
Dsp = Vsp/V;
Dtp = Vtp/V;

M = (Vq'*diag(wq)*Vq);
Pq = M\(Vq'*diag(wq));
[Vf] = pyr_basis(N,rf,sf,tf);
Ef = Vf*Pq;
Qr = Vq'*diag(wq)*Vrq;
Br = diag(nrJq.*wf);
Qs = Vq'*diag(wq)*Vsq;
Bs = diag(nsJq.*wf);
Qt = Vq'*diag(wq)*Vtq;
Bt = diag(ntJq.*wf);

norm(Qr+Qr' - Vf'*Br*Vf,'fro')
norm(Qs+Qs' - Vf'*Bs*Vf,'fro')
norm(Qt+Qt' - Vf'*Bt*Vf,'fro')

% projected div ops
Dr = M\Qr;
Ds = M\Qs;
Dt = M\Qt;

% actual nodal diff ops (exact face traces)
Dr = Vr/V;
Ds = Vs/V;
Dt = Vt/V;

% warp and make geometric terms
if 0
    a = .5;
    x = r;
    y = s;
    z = t;
    x = x + a*cos(x).*sin(y).*sin(z);
    y = y + a*sin(x).*cos(y).*sin(z);
    z = z + a*sin(x).*sin(y).*cos(z);
else
    a = .1;
    [r1 s1 t1] = pyr_nodes(1);
    V1 = pyr_basis(1,r1,s1,t1);
    VN = pyr_basis(1,r,s,t)/V1;       
    x = VN*(r1+a*(1-s1).*(1-t1));
    y = VN*(s1-a*(1-r1).*(1-t1));
    z = VN*(t1+a*(1-r1).*(1-s1));
end


if 0
    Fr = (Dr*y).*z;
    Fs = (Ds*y).*z;
    Ft = (Dt*y).*z;
    rxJ = Dt*(Fs) - Ds*(Ft);
    sxJ = Dr*(Ft) - Dt*(Fr);
    txJ = Ds*(Fr) - Dr*(Fs);
    
    Fr = (Dr*x).*z;
    Fs = (Ds*x).*z;
    Ft = (Dt*x).*z;
    ryJ = -(Dt*(Fs) - Ds*(Ft));
    syJ = -(Dr*(Ft) - Dt*(Fr));
    tyJ = -(Ds*(Fr) - Dr*(Fs));
    
    Fr = (Dr*y).*x;
    Fs = (Ds*y).*x;
    Ft = (Dt*y).*x;
    rzJ = -(Dt*(Fs) - Ds*(Ft));
    szJ = -(Dr*(Ft) - Dt*(Fr));
    tzJ = -(Ds*(Fr) - Dr*(Fs));
else
    Dr = Drp;
    Ds = Dsp;
    Dt = Dtp;
    xr = Dr*x; xs = Ds*x; xt = Dt*x;
    yr = Dr*y; ys = Ds*y; yt = Dt*y;
    zr = Dr*z; zs = Ds*z; zt = Dt*z;
    rxJ =  (ys.*zt - zs.*yt); ryJ = -(xs.*zt - zs.*xt); rzJ = (xs.*yt - ys.*xt);
    sxJ = -(yr.*zt - zr.*yt); syJ =  (xr.*zt - zr.*xt); szJ = -(xr.*yt - yr.*xt);
    txJ =  (yr.*zs - zr.*ys); tyJ = -(xr.*zs - zr.*xs); tzJ = (xr.*ys - yr.*xs);
    J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
end



nxJ = nrJ.*rxJ(Fmask) + nsJ.*sxJ(Fmask) + ntJ.*txJ(Fmask);
nyJ = nrJ.*ryJ(Fmask) + nsJ.*syJ(Fmask) + ntJ.*tyJ(Fmask);
nzJ = nrJ.*rzJ(Fmask) + nsJ.*szJ(Fmask) + ntJ.*tzJ(Fmask);

% quiver3(x(Fmask),y(Fmask),z(Fmask),nxJ,nyJ,nzJ)
hold on
% plot3(x,y,z,'o')

color_line3(ap,bp,cp,rxJ,'.'); view(3)

Nfpt = (N+1)*(N+2)/2;

% fprintf('GCL err = %g\n',norm(Qr*rxJ + Qs*sxJ + Qt*txJ,'fro'))
% norm(nxJ(1:Nfpt)-nxJ(1))


% clf
% Dr(abs(Dr)<1e-7) = 0;
% Ds(abs(Ds)<1e-7) = 0;
% Dt(abs(Dt)<1e-7) = 0;
% % [~,ids] = find(Dr(fids1,:));
% % [~,ids] = find(Dt(fids1,:));
% [~,ids] = find(Dr(fids5,:));
% % [~,ids] = find(Ds(fids5,:));
% ids = uniquetol(ids);
% plot3(r,s,t,'o')
% hold on
% plot3(r(ids),s(ids),t(ids),'r.','markersize',16)

