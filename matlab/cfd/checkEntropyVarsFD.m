clear
gamma = 1.4;

rho = 2 + rand; 
u = rand;
p = 1 + rand;

m = rho*u;
E = rho*(p/(gamma-1) + u^2/2);

% check entropy vars = d(entropy)/dU
p = @(rho,m,E) (gamma-1)*(E - .5*m^2/rho);
rhoe = @(rho,m,E) E - .5*m.^2./rho;
s = @(rho,m,E) log((gamma-1).*rhoe(rho,m,E)./(rho.^gamma));
S = @(rho,m,E) -rho.*s(rho,m,E);

delta = 1e-4;
v1 = (S(rho+delta,m,E) - S(rho-delta,m,E)) / (2*delta);
v2 = (S(rho,m+delta,E) - S(rho,m-delta,E)) / (2*delta);
v3 = (S(rho,m,E+delta) - S(rho,m,E-delta)) / (2*delta);


V1 = @(rho,m,E) (-E + rhoe(rho,m,E).*(gamma + 1 - s(rho,m,E)))./(rhoe(rho,m,E));
V2 = @(rho,m,E) (m)./(rhoe(rho,m,E));
V3 = @(rho,m,E) (-rho)./(rhoe(rho,m,E));

[v1 V1(rho,m,E)]
[v2 V2(rho,m,E)]
[v3 V3(rho,m,E)]


sV = @(V1,V2,V3) gamma - V1 + V2.^2./(2*V3);
rhoeV  = @(V1,V2,V3) ((gamma-1)./((-V3).^gamma)).^(1/(gamma-1)).*exp(-sV(V1,V2,V3)/(gamma-1));
U1 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(-V3);
U2 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(V2);
U3 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(1-V2.^2./(2*V3));

% Legendre entropy
T = @(V1,V2,V3) U1(V1,V2,V3).*V1 + U2(V1,V2,V3).*V2 + U3(V1,V2,V3).*V3 - S(U1(V1,V2,V3),U2(V1,V2,V3),U3(V1,V2,V3));

v1 = V1(rho,m,E);
v2 = V2(rho,m,E);
v3 = V3(rho,m,E);
u1 = (T(v1+delta,v2,v3) - T(v1-delta,v2,v3)) / (2*delta);
u2 = (T(v1,v2+delta,v3) - T(v1,v2-delta,v3)) / (2*delta);
u3 = (T(v1,v2,v3+delta) - T(v1,v2,v3-delta)) / (2*delta);

[u1 U1(v1,v2,v3)]
[u2 U2(v1,v2,v3)]
[u3 U3(v1,v2,v3)]
