% VXYZ = vertex positions
% r s t = values at which to evaluate geom factors

function [x,y,z,rx,sx,tx,ry,sy,ty,rz,sz,tz,J] = wedge_geom_factors(VX,VY,VZ,r,s,t)

% get vertex nodal bases 
% consistent with JF ordering
if 0
    u = [-1 1 -1 -1 1 -1]; v = [-1 -1 1 -1 -1 1]; w = [-1 -1 -1 1 1 1];
    r1 = v; s1 = w; t1 = u; % flipping coordinates for Gmsh
else
    [r1 s1 t1] = wedge_nodes(1);
end

% figure
% plot3(r1,s1,t1,'.')
% text(r1+.025,s1+.025,t1+.025,num2str((1:6)'))
% return

% r1 = [-1 1 1 -1 -1 -1]'; s1 = [-1 -1 1 1 -1 1]'; t1 = [-1 -1 -1 -1 1 1]';
V1 = wedge_basis(1,r1,s1,t1); 
invV = inv(V1);

[V1 Dr1 Ds1 Dt1] = wedge_basis(1,r,s,t); % eval map @ cubature
Interp = V1*invV; Dr = Dr1*invV; Ds = Ds1*invV; Dt = Dt1*invV;

x = Interp*VX; y = Interp*VY; z = Interp*VZ;

xr = Dr*VX; yr = Dr*VY; zr = Dr*VZ;
xs = Ds*VX; ys = Ds*VY; zs = Ds*VZ;
xt = Dt*VX; yt = Dt*VY; zt = Dt*VZ;

J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
rx =  (ys.*zt - zs.*yt)./J; ry = -(xs.*zt - zs.*xt)./J; rz = (xs.*yt - ys.*xt)./J;
sx = -(yr.*zt - zr.*yt)./J; sy =  (xr.*zt - zr.*xt)./J; sz = -(xr.*yt - yr.*xt)./J;
tx =  (yr.*zs - zr.*ys)./J; ty = -(xr.*zs - zr.*xs)./J; tz = (xr.*ys - yr.*xs)./J;

