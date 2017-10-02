function [x,y,z,rx,sx,tx,ry,sy,ty,rz,sz,tz,J] = pyr_geom_factors(VX,VY,VZ,r,s,t)

NN = 1:10;
Np = (NN+1).*(NN+2).*(2*NN+3)/6;
N = find(Np==length(VX));

[r1 s1 t1] = pyr_nodes(N);
V1 = pyr_basis(N,r1,s1,t1);
invV = inv(V1);
[V1 Dr1 Ds1 Dt1] = pyr_basis(N,r,s,t); % eval map @ cubature

if 0
    % get vertex nodal bases
    r1 = [ -1   1   1  -1  -1 ]';s1 = [ -1  -1   1   1  -1 ]';
    t1 = [ -1  -1  -1  -1   1 ]';
    
    % plot3(r1,s1,t1,'.')
    % text(r1,s1,t1,num2str((1:5)'))
    % return
    
    V1 = pyr_basis(1,r1,s1,t1);
    invV = inv(V1);
    [V1 Dr1 Ds1 Dt1] = pyr_basis(1,r,s,t); % eval map @ cubature
end

Interp = V1*invV; Dr = Dr1*invV; Ds = Ds1*invV; Dt = Dt1*invV;

x = Interp*VX; y = Interp*VY; z = Interp*VZ;

xr = Dr*VX; yr = Dr*VY; zr = Dr*VZ;
xs = Ds*VX; ys = Ds*VY; zs = Ds*VZ;
xt = Dt*VX; yt = Dt*VY; zt = Dt*VZ;

J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
rx =  (ys.*zt - zs.*yt)./J; ry = -(xs.*zt - zs.*xt)./J; rz = (xs.*yt - ys.*xt)./J;
sx = -(yr.*zt - zr.*yt)./J; sy =  (xr.*zt - zr.*xt)./J; sz = -(xr.*yt - yr.*xt)./J;
tx =  (yr.*zs - zr.*ys)./J; ty = -(xr.*zs - zr.*xs)./J; tz = (xr.*ys - yr.*xs)./J;

