[r w] = JacobiGL(0,0,N);
V = Vandermonde1D(N,r);
Vr = GradVandermonde1D(N,r);
Dr = Vr/V;
% M = inv(V*V'); 
M = diag(w); 
f = randn(N+1,1);
J = rand(N+1,1)+.5;

fb = zeros(N+1,1); fb(1) = -f(1)*J(1); fb(N+1) = f(N+1)*J(N+1);

rhs1 = diag(J)*M*Dr*f;
rhs2 = fb-diag(J)*Dr'*M*f;

rhs1-rhs2