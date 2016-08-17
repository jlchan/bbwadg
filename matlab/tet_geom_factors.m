function [rx,sx,tx,ry,sy,ty,rz,sz,tz,J] = tet_geom_factors(VX,VY,VZ,Dr,Ds,Dt)

xr = Dr*VX; yr = Dr*VY; zr = Dr*VZ;
xs = Ds*VX; ys = Ds*VY; zs = Ds*VZ;
xt = Dt*VX; yt = Dt*VY; zt = Dt*VZ;

J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
rx =  (ys.*zt - zs.*yt)./J; ry = -(xs.*zt - zs.*xt)./J; rz = (xs.*yt - ys.*xt)./J;
sx = -(yr.*zt - zr.*yt)./J; sy =  (xr.*zt - zr.*xt)./J; sz = -(xr.*yt - yr.*xt)./J;
tx =  (yr.*zs - zr.*ys)./J; ty = -(xr.*zs - zr.*xs)./J; tz = (xr.*ys - yr.*xs)./J;

