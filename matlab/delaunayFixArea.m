function tri = delaunayFixArea(vx,vy)

% delaunay w/removed small values
tri = delaunay(vx,vy);
bad = zeros(size(tri,1),1);
for e = 1:size(tri,1)
    v = [vx(tri(e,:));vy(tri(e,:))];    
    a = norm(v(:,2)-v(:,1));    b = norm(v(:,3)-v(:,2));
    c = norm(v(:,1)-v(:,3));    s = .5*(a+b+c);
    A = sqrt ( s* ( s-a)* ( s-b)* ( s-c) );
    if abs(A)<1e-5
        bad(e) = 1;
    end
end
tri(find(bad),:) = [];