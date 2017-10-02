% converts biunit to symmetric pyramid nodes 
function [r s t] = biunit_to_sym(r,s,t)

% sym verts
rs = [-1 1 1 -1 0]; ss = [-1 -1 1 1 0]; ts = [0 0 0 0 1];
rs = rs(:); ss = ss(:); ts = ts(:);
srst = [rs ss ts];

% biunit verts
rb = [-1 1 1 -1 -1]; sb = [-1 -1 1 1 -1]; tb = [-1 -1 -1 -1 1];
rb = rb(:); sb= sb(:); tb = tb(:);

% make vertex map
V1 = pyr_basis(1,rb,sb,tb); % eval map @ cubature
map = V1\srst;

% new coordinates
V = pyr_basis(1,r,s,t); % eval map @ cubature
r = V*map(:,1); s = V*map(:,2); t = V*map(:,3);

