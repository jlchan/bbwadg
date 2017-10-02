% equi nodes x equi points
function [r s t] = wedge_equi_nodes(N)

% triangle nodes
[rt st] = EquiNodes2D(N); [rt st] = xytors(rt,st); rt = rt(:); st = st(:);

% tensor product with 1D SEM nodes
%[t1D] = JacobiGL(0,0,N); 
[t1D] = linspace(-1,1,N+1);
[te re] = meshgrid(t1D,rt); [~,se] = meshgrid(t1D,st);
r = re(:); t = se(:); s = te(:);
