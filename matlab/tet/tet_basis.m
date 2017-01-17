function [V Vr Vs Vt] = tet_basis(N,r,s,t)

r = r(:); s = s(:); t = t(:);

% use codes3D
V = Vandermonde3D(N,r,s,t);
[Vr Vs Vt] = GradVandermonde3D(N,r,s,t);
