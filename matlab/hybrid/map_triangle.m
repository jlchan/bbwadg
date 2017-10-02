% map points on a reference triangle to 3D triangle plane

function [r s t] = map_triangle(VX,VY,VZ,r,s)

r1 = [-1 1 -1]'; s1 = [-1 -1 1]';
invV = inv(Vandermonde2D(1,r1,s1));

Interp = Vandermonde2D(1,r,s)*invV;

r = Interp*VX; s = Interp*VY; t = Interp*VZ;

