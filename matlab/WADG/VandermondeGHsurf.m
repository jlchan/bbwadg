% Gordon Hall blending VDM for interior nodes given surface nodes

function V = VandermondeGHsurf(N,r,s,t)

V(:,1) = -(1+r+s+t)/2;
V(:,2) = (1+r)/2;
V(:,3) = (1+s)/2;
V(:,4) = (1+t)/2;

sk = 5;

edges = [1 4; 2 4; 3 4; 1 2; 2 3; 3 1];
for e = 1:6        
    for j=0:N-2
        v1 = edges(e,1);     v2 = edges(e,2);
        t = V(:,v1)-V(:,v2);
        V(:,sk) = V(:,v1).*V(:,v2).*JacobiP(t, 1, 1, j);
        sk = sk+1;
    end
end

vv1 = [-1;-1];
vv2 = [1; -1];
vv3 = [-1; 1];
faces = [1 2 4; 2 3 4; 3 1 4; 1 2 3];
Nb = (N-3);
for f = 1:4
    v1 = faces(f,1); v2 = faces(f,2); v3 = faces(f,3);
    rr = vv1(1)*V(:,v1) + vv2(1)*V(:,v2) + vv3(1)*V(:,v3);
    ss = vv1(2)*V(:,v1) + vv2(2)*V(:,v2) + vv3(2)*V(:,v3);    
    Vf = Vandermonde2D(N-3,rr,ss);
    
    for i = 1:(Nb+1)*(Nb+2)/2 % inner face nodes
        V(:,sk) = V(:,v1).*V(:,v2).*V(:,v3).*Vf(:,i);
        sk = sk + 1;
    end
end

return;
