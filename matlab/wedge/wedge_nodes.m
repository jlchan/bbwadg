function [r s t] = wedge_nodes(N)

Nfp = (N+1)*(N+2)/2;

% get tet nodes
[x y z] = Nodes3D(N);
[r s t] = xyztorst(x,y,z);

% tet face nodes
ids = abs(s+1) < 1e-12; 
rf = r(ids); tf = t(ids);
if nnz(ids)~= Nfp
    error('Wrong # of face nodes')
end

s1D = JacobiGL(0,0,N);
% s1D = JacobiGQ(0,0,N);

r = zeros(Nfp*(N+1),1);
s = zeros(Nfp*(N+1),1);
t = zeros(Nfp*(N+1),1);
ids = 1:Nfp;
for i = 0:N
    r(ids + i*Nfp) = rf;
    s(ids + i*Nfp) = s1D(i+1)*ones(Nfp,1);
    t(ids + i*Nfp) = tf;
end
r = r(:); s = s(:); t = t(:);
% plot3(r,s,t,'.')
