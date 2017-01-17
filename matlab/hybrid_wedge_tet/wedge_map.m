function [x y z J geo] = wedge_map(r,s,t,VX,VY,VZ)

VX = VX(:); VY = VY(:); VZ = VZ(:);

%[V1 V1r V1s V1t] = bern_wedge(1,r,s,t);
[V1 V1r V1s V1t] = wedge_sem_basis(1,r,s,t);
x = V1*VX; y = V1*VY; z = V1*VZ;

xr = V1r*VX; xs = V1s*VX; xt = V1t*VX;
yr = V1r*VY; ys = V1s*VY; yt = V1t*VY;
zr = V1r*VZ; zs = V1s*VZ; zt = V1t*VZ;

J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
rx =  (ys.*zt - zs.*yt)./J; ry = -(xs.*zt - zs.*xt)./J; rz = (xs.*yt - ys.*xt)./J;
sx = -(yr.*zt - zr.*yt)./J; sy =  (xr.*zt - zr.*xt)./J; sz = -(xr.*yt - yr.*xt)./J;
tx =  (yr.*zs - zr.*ys)./J; ty = -(xr.*zs - zr.*xs)./J; tz = (xr.*ys - yr.*xs)./J;

geo.rx = rx; geo.ry = ry; geo.rz = rz;
geo.sx = sx; geo.sy = sy; geo.sz = sz;
geo.tx = tx; geo.ty = ty; geo.tz = tz;
geo.J = J;