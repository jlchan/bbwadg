N = 4;
r = JacobiGL(0,0,N);
[rq wq] = JacobiGQ(0,0,N+1);
rq = [-1;rq;1]; wq = [0;wq;0];

V = Vandermonde1D(N,r);
Vq = Vandermonde1D(N,rq)/V;
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));

rp = linspace(-1,1,250)';
Vp = Vandermonde1D(N,rp)/V;

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


dUdV11 = @(V1,V2,V3)-(V3.*exp(-(-V1+gamma+(V2.^2.*(1.0./2.0))./V3)./(gamma-1.0)).*((-V3).^(-gamma).*(gamma-1.0)).^(1.0./(gamma-1.0)))./(gamma-1.0);
dUdV12 = @(V1,V2,V3)(V2.*exp(-(-V1+gamma+(V2.^2.*(1.0./2.0))./V3)./(gamma-1.0)).*((-V3).^(-gamma).*(gamma-1.0)).^(1.0./(gamma-1.0)))./(gamma-1.0);
dUdV13 = @(V1,V2,V3)(exp((V1.*V3-V3.*gamma-V2.^2.*(1.0./2.0))./(V3.*(gamma-1.0))).*((-V3).^(-gamma).*(gamma-1.0)).^(1.0./(gamma-1.0)).*(V3.*2.0-V2.^2).*(1.0./2.0))./(V3.*(gamma-1.0));
dUdV22 = @(V1,V2,V3)-(exp((V1.*V3-V3.*gamma-V2.^2.*(1.0./2.0))./(V3.*(gamma-1.0))).*((-V3).^(-gamma).*(gamma-1.0)).^(1.0./(gamma-1.0)).*(V3-V3.*gamma+V2.^2))./(V3.*(gamma-1.0));
dUdV23 = @(V1,V2,V3)(V2.*1.0./V3.^2.*exp((V1.*V3-V3.*gamma-V2.^2.*(1.0./2.0))./(V3.*(gamma-1.0))).*((-V3).^(-gamma).*(gamma-1.0)).^(1.0./(gamma-1.0)).*(V3.*gamma.*2.0-V2.^2).*(-1.0./2.0))./(gamma-1.0);
dUdV33 = @(V1,V2,V3)(1.0./V3.^3.*exp((V1.*V3-V3.*gamma-V2.^2.*(1.0./2.0))./(V3.*(gamma-1.0))).*((-V3).^(-gamma).*(gamma-1.0)).^(1.0./(gamma-1.0)).*(V3.^2.*gamma.*4.0+V2.^4-V2.^2.*V3.*gamma.*4.0).*(-1.0./4.0))./(gamma-1.0);

% project some solution
k = pi;
rho = 2.25 + Pq*sin(k*rq);
m = Pq*sin(k*rq);
E = 2 + Pq*sin(k*rq);

rhoq = Vq*rho; 
mq = Vq*m;
Eq = Vq*E;
q1 = Pq*V1(rhoq,mq,Eq);
q2 = Pq*V2(rhoq,mq,Eq);
q3 = Pq*V3(rhoq,mq,Eq);

rhovq = U1(Vq*q1,Vq*q2,Vq*q3);
mvq = U2(Vq*q1,Vq*q2,Vq*q3);
Evq = U3(Vq*q1,Vq*q2,Vq*q3);

rhovp = U1(Vp*q1,Vp*q2,Vp*q3);
mvp = U2(Vp*q1,Vp*q2,Vp*q3);
Evp = U3(Vp*q1,Vp*q2,Vp*q3);

[q1max id1] = max(V1(rhoq,mq,Eq)); 
[q2max id2] = max(V2(rhoq,mq,Eq));
[q3max id3] = max(V3(rhoq,mq,Eq));
a = .25;
q1max = a*max(Vq*q1)+(1-a)*q1max;
q2max = a*max(Vq*q2)+(1-a)*q2max;
q3max = a*max(Vq*q3)+(1-a)*q3max;

figure(1)
hold on

plot(rp,V1(Vp*rho,Vp*m,Vp*E))
plot(rp,V2(Vp*rho,Vp*m,Vp*E))
plot(rp,V3(Vp*rho,Vp*m,Vp*E))
plot(rq,q1max*ones(size(rq)),'--.')
plot(rq,q2max*ones(size(rq)),'--.')
plot(rq,q3max*ones(size(rq)),'--.')


rhoavg = wq'*(Vq*rho)/2;
mavg = wq'*(Vq*m)/2;
Eavg = wq'*(Vq*E)/2;
t0 = 1;
t1 = t0;
t2 = .75;
t3 = .75;
rhot = rhoavg + t1*(rho-rhoavg);
mt = mavg + t2*(m-mavg);
Et = Eavg + t3*(E-Eavg);

rhoq = Vq*rhot; 
mq = Vq*mt;
Eq = Vq*Et;
q1 = Pq*V1(rhoq,mq,Eq);
q2 = Pq*V2(rhoq,mq,Eq);
q3 = Pq*V3(rhoq,mq,Eq);

% determine theta 

for i = 1:5
    a = .125;
    t1 = t1 - a*Vq*Pq*((dUdV11(Vq*q1,Vq*q2,Vq*q3).*(Vq*q1-q1max) + dUdV12(Vq*q1,Vq*q2,Vq*q3).*(Vq*q2-q2max) + dUdV13(Vq*q1,Vq*q2,Vq*q3).*(Vq*q3-q3max))./(Vq*rhot-rhoavg));
    t2 = t2 - a*Vq*Pq*((dUdV12(Vq*q1,Vq*q2,Vq*q3).*(Vq*q1-q1max) + dUdV22(Vq*q1,Vq*q2,Vq*q3).*(Vq*q2-q2max) + dUdV23(Vq*q1,Vq*q2,Vq*q3).*(Vq*q3-q3max))./(Vq*mt-mavg));
    t3 = t3 - a*Vq*Pq*((dUdV13(Vq*q1,Vq*q2,Vq*q3).*(Vq*q1-q1max) + dUdV23(Vq*q1,Vq*q2,Vq*q3).*(Vq*q2-q2max) + dUdV33(Vq*q1,Vq*q2,Vq*q3).*(Vq*q3-q3max))./(Vq*Et-Eavg));
    t1 = max(0,min(t1(id1),1));
    t2 = max(0,min(t2(id2),1));
    t3 = max(0,min(t3(id3),1));
    
    % recompute with new limited vals
    rhot = rhoavg + t1*(rho-rhoavg);
    mt = mavg + t2*(m-mavg);
    Et = Eavg + t3*(E-Eavg);
    
    q1 = Pq*V1(Vq*rhot,Vq*mt,Vq*Et);
    q2 = Pq*V2(Vq*rhot,Vq*mt,Vq*Et);
    q3 = Pq*V3(Vq*rhot,Vq*mt,Vq*Et);
    tvec = [t1,t2,t3]
end

rhoq = Vq*rhot; 
mq = Vq*mt;
Eq = Vq*Et;
q1 = Pq*V1(rhoq,mq,Eq);
q2 = Pq*V2(rhoq,mq,Eq);
q3 = Pq*V3(rhoq,mq,Eq);

plot(rp,Vp*q1,'--')
plot(rp,Vp*q2,'--')
plot(rp,Vp*q3,'--')


rhovp = U1(Vp*q1,Vp*q2,Vp*q3);
mvp = U2(Vp*q1,Vp*q2,Vp*q3);
Evp = U3(Vp*q1,Vp*q2,Vp*q3);
figure
hold on
plot(rp,rhovp)
plot(rp,Vp*rho)
% plot(rp,mvp)
% plot(rp,Evp)

