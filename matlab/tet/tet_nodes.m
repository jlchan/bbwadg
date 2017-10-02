function [r s t] = tet_nodes(N)

% get tet nodes
[x y z] = Nodes3D(N); [r s t] = xyztorst(x,y,z);
