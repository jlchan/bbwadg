% Gordon Hall blending VDM for face = 1
function V = VandermondeGH(N,r,s,t)

V(:,1) = -(1+r+s+t)/2;
V(:,2) = (1+r)/2;
V(:,3) = (1+s)/2;
V(:,4) = (1+t)/2;

sk = 5;
for i = 0:N-2 % edge basis
    V(:,sk) = V(:,1).*V(:,2).*JacobiP(V(:,1)-V(:,2),0,0,i); sk = sk + 1;
    V(:,sk) = V(:,2).*V(:,3).*JacobiP(V(:,2)-V(:,3),0,0,i); sk = sk + 1;
    V(:,sk) = V(:,3).*V(:,1).*JacobiP(V(:,3)-V(:,1),0,0,i); sk = sk + 1;
end

Nb = (N-3);
Vface = Vandermonde2D(Nb,r,s); % lower degree polynomials
for i = 1:(Nb+1)*(Nb+2)/2 % face nodes
    V(:,sk) = V(:,1).*V(:,2).*V(:,3).*Vface(:,i);
    sk = sk + 1;
end