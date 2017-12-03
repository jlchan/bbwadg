clear
N = 5;

global gamma
gamma = 1.4;

rhoe = @(rho,rhou,rhov,E) E - .5*(rhou.^2+rhov.^2)./rho;
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

Nq = 4*N;
[rq sq wq] = Cubature2D(Nq);
[r s] = Nodes2D(N); [r s] = xytors(r,s);
V = Vandermonde2D(N,r,s);
Vq = Vandermonde2D(N,rq,sq)/V;
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));

rho = Pq*(2 + sin(pi*rq).*sin(pi*sq));
rhou = .1*randn(size(r));
rhov = .1*randn(size(r));
E = Pq*(1.25 + cos(pi*rq).*cos(pi*sq));

rhoavg = V*diag([1;zeros(length(r)-1,1)])*(V\rho);
rhouavg = V*diag([1;zeros(length(r)-1,1)])*(V\rhou);
rhovavg = V*diag([1;zeros(length(r)-1,1)])*(V\rhov);
Eavg = V*diag([1;zeros(length(r)-1,1)])*(V\E);

trho = 1; rho = rhoavg + trho*(rho-rhoavg);
trhou = 1; rhou = rhouavg + trhou*(rhou-rhouavg);
trhov = 1; rhov = rhovavg + trhov*(rhov-rhovavg);
tE = .1; E = Eavg + tE*(E-Eavg);


q1 = Vq*Pq*V1(Vq*rho,Vq*rhou,Vq*rhov,Vq*E);
q2 = Vq*Pq*V2(Vq*rho,Vq*rhou,Vq*rhov,Vq*E);
q3 = Vq*Pq*V3(Vq*rho,Vq*rhou,Vq*rhov,Vq*E);
q4 = Vq*Pq*V4(Vq*rho,Vq*rhou,Vq*rhov,Vq*E);

rhoq = U1(q1,q2,q3,q4);
rhouq = U2(q1,q2,q3,q4);
rhovq = U3(q1,q2,q3,q4);
Eq = U4(q1,q2,q3,q4);

vv = Vq*rho;
color_line3(rq,sq,vv,vv,'.')
hold on

vv = rhoq;
color_line3(rq,sq,vv,vv,'o')
view(3)
