% test sem-lumped BB lift

N = 3;

[r w] = JacobiGQ(0,0,N);
[rSEM] = JacobiGL(0,0,N); V = Vandermonde1D(N,rSEM); wSEM = sum(inv(V*V'),2);
r = rSEM; w = wSEM; % use SEM quadrature

V = bern_basis_1D(N,r);
% V = eye(N+1);

J = rand*((r+1)/2) + rand;
% J = 1;
M = V'*diag(w.*J)*V;
Mf = [1;zeros(N,1)];

LIFT = M\Mf;
% min(abs(LIFT(LIFT>1e-8))).*w(1).*J(1)
% LIFT = LIFT/min(abs(LIFT(LIFT>1e-8)));
LIFT*(w(1)*J(1))
return

req = linspace(-1,1,N+1);
hold on;plot(req,LIFT,'o-')
rp = linspace(-1,1,250);
Vp = bern_basis_1D(N,rp);
plot(rp,Vp*LIFT)

rGR = JacobiGL(0,0,N);
plot(rp,0*rp,'k--')
plot(rGR,rGR*0,'s','markersize',8)