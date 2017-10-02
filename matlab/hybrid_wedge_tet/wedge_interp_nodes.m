% WB nodes x SEM points
function [r s t] = wedge_interp_nodes(N)

% triangle nodes
[rt st] = Nodes2D(N); [rt st] = xytors(rt,st); rt = rt(:); st = st(:);

% tensor product with 1D SEM nodes
[t1D] = JacobiGL(0,0,N); 
[te re] = meshgrid(t1D,rt); [~,se] = meshgrid(t1D,st);
r = re(:); s = se(:); t = te(:);
