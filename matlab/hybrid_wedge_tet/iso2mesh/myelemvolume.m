function vol = myelemvolume(node,elem)

r = node(:,1); s = node(:,2); t = node(:,3);

vol = zeros(size(elem,1),1);
for k=1:size(elem,1)
    vids = elem(k,1:4);
    V = [r(vids(:)) s(vids(:)) t(vids(:)) ones(4,1)];        
    vol(k) = abs(det(V)/6);
end

