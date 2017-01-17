% map points on a reference triangle to 3D triangle plane

function [r s t] = map_quad(VX,VY,VZ,r,s)

r1 = [-1 1 1 -1]'; s1 = [-1 -1 1 1]';
invV = inv(quad_basis(1,r1,s1));

Interp = quad_basis(1,r,s)*invV;

r = Interp*VX; s = Interp*VY; t = Interp*VZ;

function V = quad_basis(N, r, s)

sk = 1;
V = zeros(length(r), (N+1)*(N+1));
for i=0:N
    for j=0:N
        V(:,sk) = JacobiP(r, 0, 0, i).*JacobiP(s, 0, 0, j);
        sk = sk+1;
    end
end