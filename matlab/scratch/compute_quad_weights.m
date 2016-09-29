clear

N = 9;

[r s w] = Cubature2D(2*N+1);

V = Vandermonde2D(N,r,s);
Np = size(V,2);

length(r)-Np
return
[~,~,E] = qr(V',0);

E = E(1:Np);
plot(r,s,'o')
hold on

r = r(E); s = s(E);
plot(r,s,'x')

V = Vandermonde2D(N,r,s);
cond(V)

% [r s] = Nodes2D(N); [r s] = xytors(r,s);
% V = Vandermonde2D(N,r,s);
% cond(V)

Vq = Vandermonde2D(2*N,r,s);
b = zeros(size(Vq',1),1); b(1) = sqrt(2);
w2 = (Vq')\b;
norm(Vq'*w2-b)

