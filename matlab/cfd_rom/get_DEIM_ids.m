function p = get_DEIM_ids(V)

VDEIM(:,1) = V(:,1);
[~,id] = max(abs(V(:,1)));
p = id;
for j = 2:size(V,2)
    r = V(:,j)-VDEIM*(VDEIM(p,:)\V(p,j));
    [~,id] = max(abs(r));
    p(j) = id;
    VDEIM = [VDEIM r];    
end
p = sort(p);
% rq = rqref(p);