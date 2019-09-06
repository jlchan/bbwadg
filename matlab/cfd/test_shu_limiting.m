N = 7;
r = JacobiGL(0,0,N);
[rq wq] = JacobiGQ(0,0,N+2);
rp = linspace(-1,1,1000)';

V = Vandermonde1D(N,r);
Vq = Vandermonde1D(N,rq)/V;
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));
Vp = Vandermonde1D(N,rp)/V;

%% swe

h = 2 + rand(N+1,1);
hv = 5*randn(N+1,1);

g = 1;

hq = Vq*h;
hvq = Vq*hv;
vq = hvq./hq;

% project entropy variables
q1 = Vq*Pq*(g*hq-vq.^2/2);
q2 = Vq*Pq*(vq);

% redefine flux variables
hq1 = (q1 + q2.^2/2)/g; % 1/g*[(g*hq-P(hv^2/h^2)/2) + P(vq)^2/2]
hq = Vq*h + .5*((Vq*Pq*vq).^2 - Vq*Pq*(vq.^2));
vq = q2;

plot(rp,Vp*h)
hold on
plot(rq,hq,'o--')

wq'*(Vq*h)-wq'*hq % avg(h) >= avg(he) = avg(h) + .5*((P*v)^2 - P*(v^2))

%% Euler

gamma = 1.4;
rhoe = @(rho,m,E) E - .5*m.^2./rho;
s = @(rho,m,E) log((gamma-1).*rhoe(rho,m,E)./(rho.^gamma));

V1 = @(rho,m,E) (-E + rhoe(rho,m,E).*(gamma + 1 - s(rho,m,E)))./(rhoe(rho,m,E));
V2 = @(rho,m,E) (m)./(rhoe(rho,m,E));
V3 = @(rho,m,E) (-rho)./(rhoe(rho,m,E));

sV = @(V1,V2,V3) gamma - V1 + V2.^2./(2*V3);
rhoeV  = @(V1,V2,V3) ((gamma-1)./((-V3).^gamma)).^(1/(gamma-1)).*exp(-sV(V1,V2,V3)/(gamma-1));
U1 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(-V3);
U2 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(V2);
U3 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(1-V2.^2./(2*V3));

rho = 2 + rand(N+1,1);
rhou = .5*randn(N+1,1);
E = rho.^gamma;

q1 = V1(Vq*rho,Vq*rhou,Vq*E);
q2 = V2(Vq*rho,Vq*rhou,Vq*E);
q3 = V3(Vq*rho,Vq*rhou,Vq*E);

u1 = U1(Vq*Pq*q1,Vq*Pq*q2,Vq*Pq*q3);
u2 = U2(Vq*Pq*q1,Vq*Pq*q2,Vq*Pq*q3);
u3 = U3(Vq*Pq*q1,Vq*Pq*q2,Vq*Pq*q3);

wq'*(Vq*rho) - wq'*u1
wq'*(Vq*rhou) - wq'*u2
wq'*(Vq*E) - wq'*u3

