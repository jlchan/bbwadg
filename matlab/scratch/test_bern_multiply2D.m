function test_bern_multiply2D

N = 2;
Np = (N+1)*(N+2)/2;
u = (1:Np)';
v = zeros(Np,1);
v(1) = 1;
% v = u;

[r s] = Nodes2D(2*N); [r s] = xytors(r,s);
[rq sq w] = Cubature2D(4*N);
V = bern_basis_tri(N,r,s);

V2N = bern_basis_tri(2*N,r,s);
uv = V2N\((V*u).*(V*v));

norm(V2N*uv-(V*u).*(V*v))

cN = bern_coeffs(N);
c2N = bern_coeffs(2*N);
Vs = bern_basis_tri(N,r,s)*diag(1./cN);
uc = u.*cN;
vc = v.*cN;
norm(V*u - Vs*uc)

keyboard


function c = bern_coeffs(N)
Np = (N+1)*(N+2)/2;

c = zeros(Np,1);
sk = 1;
for i = 0:N
    for j = 0:N-i
        k = N-i-j;
        c(sk) = factorial(N)/(factorial(i)*factorial(j)*factorial(k));
        sk = sk + 1;
    end
end
