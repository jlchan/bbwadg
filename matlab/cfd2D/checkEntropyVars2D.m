gamma = 1.4;
rho = 2;
u = 1;
v = .5;
rhou = rho.*u;
rhov = rho.*v;
E = 3;

rhoe = @(rho,rhou,rhov,E) E - (rhou.^2+rhov.^2)./(2*rho);
s = @(rho,rhou,rhov,E) log((gamma-1)*rhoe(rho,rhou,rhov,E)./(rho.^gamma));
V1 = @(rho,rhou,rhov,E) (-E + rhoe(rho,rhou,rhov,E).*(gamma + 1 - s(rho,rhou,rhov,E)))./(rhoe(rho,rhou,rhov,E));
V2 = @(rho,rhou,rhov,E) rhou./(rhoe(rho,rhou,rhov,E));
V3 = @(rho,rhou,rhov,E) rhov./(rhoe(rho,rhou,rhov,E));
V4 = @(rho,rhou,rhov,E) (-rho)./(rhoe(rho,rhou,rhov,E));

sV = @(V1,V2,V3,V4) gamma - V1 + (V2.^2+V3.^2)./(2*V4);
rhoeV  = @(V1,V2,V3,V4) ((gamma-1)./((-V4).^gamma)).^(1/(gamma-1)).*exp(-sV(V1,V2,V3,V4)/(gamma-1));
U1 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(-V4);
U2 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(V2);
U3 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(V3);
U4 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(1-(V2.^2+V3.^2)./(2*V4));

v1 = V1(rho,rhou,rhov,E);
v2 = V2(rho,rhou,rhov,E);
v3 = V3(rho,rhou,rhov,E);
v4 = V4(rho,rhou,rhov,E);
norm(rho-U1(v1,v2,v3,v4))
norm(rhou-U2(v1,v2,v3,v4))
norm(rhov-U3(v1,v2,v3,v4))
norm(E-U4(v1,v2,v3,v4))



return
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
% v1*f1 + v2*f2 + V4*f3

% v1
% (gamma-s(rho,rhou,rhov,E))./(gamma-1) - beta(rho,u,E).*u.^2