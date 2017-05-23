% function [x y z geofacs mesh] = Elbow3D(r,s,t)

% h = 1; w = 1; r = .5; R = 1; t = r/2;
h = 2; w = 2; rad = 0.5; R = 2; T = 0.2; 

[N,N_2,N_3,n_1,n_2,n_3,Xi_1,Xi_2,Xi_3,P,W] = NURBS_Pipe_Elbow_Volume(h,w,rad,R,T);

d = 3;
[K,C_operators,IEN] = Extract_Basis(N,N_2,N_3,n_1,n_2,n_3,Xi_1,Xi_2,Xi_3);

[P_b,w_b] = Extract_Geometry(d,K,C_operators,IEN,P,W);

% Bezier_Plotter(N,N_2,N_3,K,P_b,w_b);
% return

r1D = linspace(-1,1,N+1)';
r1D = linspace(-1,1,5)';
[r s t] = meshgrid(r1D); r = r(:); s = s(:); t = t(:);

Nfp = length(r1D)^2;
Nfaces = 6;

Fmask = zeros(Nfp,Nfaces);
Fmask(:,1) = find(abs(r+1) <1e-8);
Fmask(:,2) = find(abs(r-1) <1e-8);
Fmask(:,3) = find(abs(s+1) <1e-8);
Fmask(:,4) = find(abs(s-1) <1e-8);
Fmask(:,5) = find(abs(t+1) <1e-8);
Fmask(:,6) = find(abs(t-1) <1e-8);

Np = (N+1)^3;
cx = reshape(P_b(:,1,:),Np,K);
cy = reshape(P_b(:,2,:),Np,K);
cz = reshape(P_b(:,3,:),Np,K);
cw = reshape(w_b,Np,K);
% cw = cw.^0; % restore polynomial approx


if 0
    [V1D Vr1D] = bern_basis_1D(N,r1D);
    V = kron(kron(V1D,V1D),V1D);
    Vr = kron(kron(V1D,Vr1D),V1D);
    Vs = kron(kron(V1D,V1D),Vr1D);
    Vt = kron(kron(Vr1D,V1D),V1D);
else
    [V1Da Vr1Da] = bern_basis_1D(N,r);
    [V1Db Vr1Db] = bern_basis_1D(N,s);
    [V1Dc Vr1Dc] = bern_basis_1D(N,t);
    a = 2;
    Vr1Da = Vr1Da*a;
    Vr1Db = Vr1Db*a;
    Vr1Dc = Vr1Dc*a;
    V = zeros(length(r),(N+1)^3);
    Vr = zeros(length(r),(N+1)^3);
    Vs = zeros(length(r),(N+1)^3);
    Vt = zeros(length(r),(N+1)^3);
    sk = 1;
    for k = 1:N+1
        for j = 1:N+1
            for i = 1:N+1
                V(:,sk) = V1Da(:,i).*V1Db(:,j).*V1Dc(:,k);
                Vr(:,sk) = Vr1Da(:,i).*V1Db(:,j).*V1Dc(:,k);
                Vs(:,sk) = V1Da(:,i).*Vr1Db(:,j).*V1Dc(:,k);
                Vt(:,sk) = V1Da(:,i).*V1Db(:,j).*Vr1Dc(:,k);
                sk = sk + 1;
            end            
        end
    end    
end

% D1D = V1D\Vr1D; D1D(abs(D1D)<1e-8) = 0; I = eye(N+1);
% Dr = kron(kron(I,D1D),I);
% Dt = kron(kron(D1D,I),I);
% Ds = kron(kron(I,I),D1D);

w = V*cw;
x = (V*(cx.*cw))./w; 
y = (V*(cy.*cw))./w;
z = (V*(cz.*cw))./w;

wr = Vr*cw; 
ws = Vs*cw; 
wt = Vt*cw; 
xr = (Vr*(cx.*cw))./w - wr.*x./w; 
xs = (Vs*(cx.*cw))./w - ws.*x./w;
xt = (Vt*(cx.*cw))./w - wt.*x./w;

yr = (Vr*(cy.*cw))./w - wr.*y./w; 
ys = (Vs*(cy.*cw))./w - ws.*y./w;
yt = (Vt*(cy.*cw))./w - wt.*y./w;

zr = (Vr*(cz.*cw))./w - wr.*z./w; 
zs = (Vs*(cz.*cw))./w - ws.*z./w;
zt = (Vt*(cz.*cw))./w - wt.*z./w;

J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
rx =  (ys.*zt - zs.*yt)./J; ry = -(xs.*zt - zs.*xt)./J; rz = (xs.*yt - ys.*xt)./J;
sx = -(yr.*zt - zr.*yt)./J; sy =  (xr.*zt - zr.*xt)./J; sz = -(xr.*yt - yr.*xt)./J;
tx =  (yr.*zs - zr.*ys)./J; ty = -(xr.*zs - zr.*xs)./J; tz = (xr.*ys - yr.*xs)./J;
rxJ =  (ys.*zt - zs.*yt); ryJ = -(xs.*zt - zs.*xt); rzJ = (xs.*yt - ys.*xt);
sxJ = -(yr.*zt - zr.*yt); syJ =  (xr.*zt - zr.*xt); szJ = -(xr.*yt - yr.*xt);
txJ =  (yr.*zs - zr.*ys); tyJ = -(xr.*zs - zr.*xs); tzJ = (xr.*ys - yr.*xs);

geofacs.rxJ = rxJ; geofacs.sxJ = sxJ; geofacs.txJ = txJ;
geofacs.rxJ = rxJ; geofacs.sxJ = sxJ; geofacs.txJ = txJ;
geofacs.rxJ = rxJ; geofacs.sxJ = sxJ; geofacs.txJ = txJ;
geofacs.J = J;
% color_line3(x,y,z,J,'.')
% axis equal
% colorbar

xf = x(Fmask(:),:);
yf = y(Fmask(:),:);
zf = z(Fmask(:),:);

nr = zeros(length(Fmask(:)),1);
ns = zeros(length(Fmask(:)),1);
nt = zeros(length(Fmask(:)),1);

Nfp = size(Fmask,1);
fids = @(f) (1:Nfp) + (f-1)*Nfp;
f = 1; nr(fids(f)) = -1; 
f = 2; nr(fids(f)) = 1;  
f = 3; ns(fids(f)) = -1; 
f = 4; ns(fids(f)) = 1;  
f = 5; nt(fids(f)) = -1; 
f = 6; nt(fids(f)) = 1;

rf = r(Fmask(:),:); sf = s(Fmask(:),:); tf = t(Fmask(:),:);

% plot3(xf,yf,zf,'o');axis equal; return
% plot3(rf,sf,tf,'o'); hold on; quiver3(rf,sf,tf,nr,ns,nt,2); return

% phys normals
nr = repmat(nr,1,K); ns = repmat(ns,1,K); nt = repmat(nt,1,K); 

rxf = rxJ(Fmask(:),:); sxf = sxJ(Fmask(:),:); txf = txJ(Fmask(:),:);
ryf = ryJ(Fmask(:),:); syf = syJ(Fmask(:),:); tyf = tyJ(Fmask(:),:);
rzf = rzJ(Fmask(:),:); szf = szJ(Fmask(:),:); tzf = tzJ(Fmask(:),:);
% rxf = rx(Fmask(:),:); sxf = sx(Fmask(:),:); txf = tx(Fmask(:),:);
% ryf = ry(Fmask(:),:); syf = sy(Fmask(:),:); tyf = ty(Fmask(:),:);
% rzf = rz(Fmask(:),:); szf = sz(Fmask(:),:); tzf = tz(Fmask(:),:);

nxJ = rxf.*nr + sxf.*ns + txf.*nt;
nyJ = ryf.*nr + syf.*ns + tyf.*nt;
nzJ = rzf.*nr + szf.*ns + tzf.*nt;

sJ = sqrt(nxJ.^2 + nyJ.^2 + nzJ.^2);
nx = nxJ./sJ;
ny = nyJ./sJ;
nz = nzJ./sJ;

% plot3(xf,yf,zf,'o');return

e = 1:K; 
fids = 1:length(Fmask(:));
% e = 1; 
% f = 4; fids = (1:Nfp) + (f-1)*Nfp;

% vv = cy(:,e).*cw(:,e); h = color_line3(cx(:,e),cy(:,e),cz(:,e),vv,'.');  set(h,'markersize',32); colorbar; return


plot3(xf(fids,e),yf(fids,e),zf(fids,e),'o');hold on
% h = color_line3(xf(fids,e),yf(fids,e),zf(fids,e),w(fids,e),'.'); set(h,'markersize',32); colorbar
quiver3(xf(fids,e),yf(fids,e),zf(fids,e),nx(fids,e),ny(fids,e),nz(fids,e))
axis equal
view(90,0)

return

re1D = [-1; 1]; % vertices
[re se te] = meshgrid(re1D); 
re = re(:); se = se(:); te = te(:);
[V1D Vr1D] = bern_basis_1D(N,re1D);

V = kron(kron(V1D,V1D),V1D);

w = V*cw;
vx = (V*(cx.*cw))./w; 
vy = (V*(cy.*cw))./w;
vz = (V*(cz.*cw))./w;


% Vmask = [7 25 27 9 1 19 21 3]'; % vertices 
% VX = cx(Vmask,:);
% VY = cy(Vmask,:);
% VZ = cz(Vmask,:);

plot3(vx,vy,vz,'o')
axis equal

[VXYZ idr idc] = uniquetol([vx(:) vy(:) vz(:)],1e-7,'ByRows',true);
VX = VXYZ(:,1); VY = VXYZ(:,2); VZ = VXYZ(:,3);
EToV = reshape(idc,8,12)';

mesh.VX = VX;
mesh.VY = VY;
mesh.VZ = VZ;
mesh.EToV = EToV;
% for e = 1:K
%     ev = EToV(e,:);
%     clf
%     plot3(VX(ev),VY(ev),VZ(ev),'o')    
%     axis equal
%     pause
% end



[V1D Vr1D] = bern_basis_1D(N,JacobiGL(0,0,N));
V = kron(kron(V1D,V1D),V1D);
invV = kron(kron(inv(V1D),inv(V1D)),inv(V1D));
x = (V*(cx.*cw))./(V*cw);

% interpolate to spline basis
cx = invV*(x);
xr = Vr*cx;
xs = Vs*cx;
xt = Vt*cx;

% norm(rx.*xr + sx.*xs + tx.*xt - 1,'fro')
 
% clf
% color_line3(x,y,z,x,'.')

