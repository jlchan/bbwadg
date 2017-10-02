function [rx,sx,tx,ry,sy,ty,rz,sz,tz,J] = geofacs(x,y,z)

[r s t] = Nodes3D(1); [r s t] = xyztorst(r,s,t);
[V] = Vandermonde3D(1,r,s,t);
[Vr Vs Vt] = GradVandermonde3D(1,r,s,t);
Dr = Vr/V; Ds = Vs/V; Dt = Vt/V;

% calculate geometric factors
xr = Dr*x; xs = Ds*x; xt = Dt*x;
yr = Dr*y; ys = Ds*y; yt = Dt*y;
zr = Dr*z; zs = Ds*z; zt = Dt*z;

J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
rx =  (ys.*zt - zs.*yt)./J; ry = -(xs.*zt - zs.*xt)./J; rz = (xs.*yt - ys.*xt)./J;
sx = -(yr.*zt - zr.*yt)./J; sy =  (xr.*zt - zr.*xt)./J; sz = -(xr.*yt - yr.*xt)./J;
tx =  (yr.*zs - zr.*ys)./J; ty = -(xr.*zs - zr.*xs)./J; tz = (xr.*ys - yr.*xs)./J;
