% clear all; clc;

h = 1; w = 1; r = .5; R = 1; t = r/2;
h = 2; w = 2; r = 0.5; R = 2; t = 0.2; 

[p_1,p_2,p_3,n_1,n_2,n_3,Xi_1,Xi_2,Xi_3,P,W] = NURBS_Pipe_Elbow_Volume(h,w,r,R,t);

d = 3;
[n_el,C_operators,IEN] = Extract_Basis(p_1,p_2,p_3,n_1,n_2,n_3,Xi_1,Xi_2,Xi_3);

[P_b,w_b] = Extract_Geometry(d,n_el,C_operators,IEN,P,W);

% Bezier_Plotter(p_1,p_2,p_3,n_el,P_b,w_b);
% return

Np = (p_1+1)^3;
cx = reshape(P_b(:,1,:),Np,n_el);
cy = reshape(P_b(:,2,:),Np,n_el);
cz = reshape(P_b(:,3,:),Np,n_el);
cw = reshape(w_b,Np,n_el);
% cw = cw.^0;

xp = linspace(-1,1,8);
[V1D Vr1D] = bern_basis_1D(p_1,xp);
V = kron(kron(V1D,V1D),V1D);
Vr = kron(kron(V1D,Vr1D),V1D);
Vs = kron(kron(Vr1D,V1D),V1D);
Vt = kron(kron(V1D,V1D),Vr1D);

% D1D = V1D\Vr1D; D1D(abs(D1D)<1e-8) = 0; I = eye(p_1+1);
% Dr = kron(kron(I,D1D),I);
% Dt = kron(kron(D1D,I),I);
% Ds = kron(kron(I,I),D1D);

r1D = linspace(-1,1,p_1+1)';
[r s t] = meshgrid(r1D); r = r(:); s = s(:); t = t(:);
w = V*cw;
x = (V*(cx.*cw))./w;
y = (V*(cy.*cw))./w;
z = (V*(cz.*cw))./w;
% plot3(x,y,z,'o')
ids = abs(z-1)<1e-8;
% plot3(x(ids),y(ids),z(ids),'o')
axis equal

% return
wr = Vr*cw;
ws = Vs*cw;
wt = Vt*cw;
xr = (Vr*(cx.*cw))./w - wr.*x./w.^2; 
xs = (Vs*(cx.*cw))./w - ws.*x./w.^2;
xt = (Vt*(cx.*cw))./w - wt.*x./w.^2;

yr = (Vr*(cy.*cw))./w - wr.*y./w.^2; 
ys = (Vs*(cy.*cw))./w - ws.*y./w.^2;
yt = (Vt*(cy.*cw))./w - wt.*y./w.^2;

zr = (Vr*(cz.*cw))./w - wr.*z./w.^2; 
zs = (Vs*(cz.*cw))./w - ws.*z./w.^2;
zt = (Vt*(cz.*cw))./w - wt.*z./w.^2;

J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
rx =  (ys.*zt - zs.*yt)./J; ry = -(xs.*zt - zs.*xt)./J; rz = (xs.*yt - ys.*xt)./J;
sx = -(yr.*zt - zr.*yt)./J; sy =  (xr.*zt - zr.*xt)./J; sz = -(xr.*yt - yr.*xt)./J;
tx =  (yr.*zs - zr.*ys)./J; ty = -(xr.*zs - zr.*xs)./J; tz = (xr.*ys - yr.*xs)./J;

color_line3(x,y,z,J,'.')
axis equal
colorbar

[V1D Vr1D] = bern_basis_1D(p_1,JacobiGL(0,0,p_1));
V = kron(kron(V1D,V1D),V1D);
invV = kron(kron(inv(V1D),inv(V1D)),inv(V1D));
x = (V*(cx.*cw))./(V*cw);

% interpolate to spline basis
cx = invV*(x);
xr = Vr*cx;
xs = Vs*cx;
xt = Vt*cx;

% norm(rx.*xr + sx.*xs + tx.*xt - 1,'fro')
rxJ =  (ys.*zt - zs.*yt); ryJ = -(xs.*zt - zs.*xt); rzJ = (xs.*yt - ys.*xt);
sxJ = -(yr.*zt - zr.*yt); syJ =  (xr.*zt - zr.*xt); szJ = -(xr.*yt - yr.*xt);
txJ =  (yr.*zs - zr.*ys); tyJ = -(xr.*zs - zr.*xs); tzJ = (xr.*ys - yr.*xs);
norm(rxJ.*xr + sxJ.*xs + txJ.*xt - J,'fro')
% clf
% color_line3(x,y,z,x,'.')

