gamma = 1.4;
rho = 2;
m = 1;
u = m./rho;
E = 3;

rhoe = @(rho,m,E) E - m.^2./(2*rho);
s = @(rho,m,E) log((gamma-1)*rhoe(rho,m,E)./(rho.^gamma));

V1 = @(rho,m,E) (-E + rhoe(rho,m,E).*(gamma + 1 - s(rho,m,E)))./(rhoe(rho,m,E));
V2 = @(rho,m,E) (m)./(rhoe(rho,m,E));
V3 = @(rho,m,E) (-rho)./(rhoe(rho,m,E));

sV = @(V1,V2,V3) gamma - V1 + V2.^2./(2*V3);
rhoeV  = @(V1,V2,V3) ((gamma-1)./((-V3).^gamma)).^(1/(gamma-1)).*exp(-sV(V1,V2,V3)/(gamma-1));
U1 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(-V3);
U2 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(V2);
U3 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(1-V2.^2./(2*V3));

v1 = V1(rho,m,E);
v2 = V2(rho,m,E);
v3 = V3(rho,m,E);
U1(v1,v2,v3)
U2(v1,v2,v3)
U3(v1,v2,v3)

% check chandreshekar flux
avg = @(uL,uR) .5*(uL+uR);
pfun = @(rho,u,E) (gamma-1)*(E-.5*rho.*u.^2);
beta = @(rho,u,E) rho./(2*pfun(rho,u,E));
F1 = @(rhoL,rhoR,uL,uR,EL,ER) logmean(rhoL,rhoR).*avg(uL,uR);
F2 = @(rhoL,rhoR,uL,uR,EL,ER) avg(rhoL,rhoR)./(2*avg(beta(rhoL,uL,EL),beta(rhoR,uR,ER))) + avg(uL,uR).*F1(rhoL,rhoR,uL,uR,EL,ER);
F3 = @(rhoL,rhoR,uL,uR,EL,ER) F1(rhoL,rhoR,uL,uR,EL,ER)...
    .*(1./(2*(gamma-1).*logmean(beta(rhoL,uL,EL),beta(rhoR,uR,ER))) - .5*avg(uL.^2,uR.^2)) ...
    + avg(uL,uR).*F2(rhoL,rhoR,uL,uR,EL,ER);

% f1 = F1(
% v1*f1 + v2*f2 + v3*f3

v1
(gamma-s(rho,m,E))./(gamma-1) - beta(rho,u,E).*u.^2