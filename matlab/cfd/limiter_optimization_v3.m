clear
N = 7;
r = JacobiGL(0,0,N);
[rq wq] = JacobiGQ(0,0,N+1);
rq = [-1;rq;1]; wq = [0;wq;0];
Nq = length(rq);

V = Vandermonde1D(N,r);
Vq = Vandermonde1D(N,rq)/V;
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));
VqPq = Vq*Pq;
rp = linspace(-1,1,250)';
Vp = Vandermonde1D(N,rp)/V;

% define functions
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

U1v = @(V,rho,m,E) U1(V*Pq*V1(rho,m,E),V*Pq*V2(rho,m,E),V*Pq*V3(rho,m,E));
U2v = @(V,rho,m,E) U2(V*Pq*V1(rho,m,E),V*Pq*V2(rho,m,E),V*Pq*V3(rho,m,E));
U3v = @(V,rho,m,E) U3(V*Pq*V1(rho,m,E),V*Pq*V2(rho,m,E),V*Pq*V3(rho,m,E));

% define vars
a = 1;
b = 1/2;
k = pi;
rho = 1 - b*.25*Pq*exp(-10*rq.^2);
m = 1 + b*Pq*sin(a+k*rq);
E = 2 + b*Pq*sin(a+k*rq);

rhoq = Vq*rho;
mq = Vq*m;
Eq = Vq*E;
rhoavg = .5*wq'*rhoq;
mavg = .5*wq'*mq;
Eavg = .5*wq'*Eq;

v1 = V1(rhoq,mq,Eq);
v2 = V2(rhoq,mq,Eq);
v3 = V3(rhoq,mq,Eq);
pv1 = Vq*Pq*v1;
pv2 = Vq*Pq*v2;
pv3 = Vq*Pq*v3;

rhovq = U1v(Vq,rhoq,mq,Eq);
mvq = U2v(Vq,rhoq,mq,Eq);
Evq = U3v(Vq,rhoq,mq,Eq);

lim = @(u,theta) .5*wq'*u + theta*(u - .5*wq'*u);
F = V*(diag([ones(N,1);.1])/V);

% minimize .5*(S(uv)-S(u))^2 + (t-1).^2
t = [1 1 1];
Vf1 = Vandermonde1D(N,1)/V;

S = @(rho,m,E) -rho.*s(rho,m,E);
% conjugate entropies 
TV = @(V1,V2,V3) U1(V1,V2,V3).*V1 + U2(V1,V2,V3).*V2 + U3(V1,V2,V3).*V3 - S(U1(V1,V2,V3),U2(V1,V2,V3),U3(V1,V2,V3));
TU = @(U1,U2,U3) U1.*V1(U1,U2,U3) + U2.*V2(U1,U2,U3) + U3.*V3(U1,U2,U3) - S(U1,U2,U3);
TTV = @(V1,V2,V3) U1(V1,V2,V3).*V1 + U2(V1,V2,V3).*V2 + U3(V1,V2,V3).*V3;
TTU = @(U1,U2,U3) U1.*V1(U1,U2,U3) + U2.*V2(U1,U2,U3) + U3.*V3(U1,U2,U3);

V30 = V3(Vf1*rho,Vf1*m,Vf1*E);
J = @(t) (Vf1*Pq*V3(lim(rhoq,t(1)),lim(mq,t(2)),lim(Eq,t(3))) - V30)^2;

% T0 = max(TU(Vq*rho,Vq*m,Vq*E));
% J = @(t) TV(Vf1*Pq*V1(lim(rhoq,t(1)),lim(mq,t(2)),lim(Eq,t(3))),...
%     Vf1*Pq*V2(lim(rhoq,t(1)),lim(mq,t(2)),lim(Eq,t(3))),...
%     Vf1*Pq*V3(lim(rhoq,t(1)),lim(mq,t(2)),lim(Eq,t(3)))- T0)^2;

% T0 = max(TTU(Vq*rho,Vq*m,Vq*E));
% J = @(t) TTV(Vf1*Pq*V1(lim(rhoq,t(1)),lim(mq,t(2)),lim(Eq,t(3))),...
%     Vf1*Pq*V2(lim(rhoq,t(1)),lim(mq,t(2)),lim(Eq,t(3))),...
%     Vf1*Pq*V3(lim(rhoq,t(1)),lim(mq,t(2)),lim(Eq,t(3)))- T0)^2;

% [t,fval,flag,out] = fminsearch(@(t) J(t) + norm(t-1)^2,t);
% t = min(t,1);
% t = [1 1 0];
% t = [1 .9 1];
% t = [.5 1 1];

hold on
plot(rp,V3(Vp*rho,Vp*m,Vp*E),'b-')
plot(rp,Vp*Pq*V3(rhoq,mq,Eq),'r-')
plot(rp,V3(Vp*Pq*lim(rhoq,t(1)),Vp*Pq*lim(mq,t(2)),Vp*Pq*lim(Eq,t(3))),'k-')
plot(rp,Vp*Pq*V3(lim(rhoq,t(1)),lim(mq,t(2)),lim(Eq,t(3))),'k--')
title(sprintf('t = [%1.2f %1.2f %1.2f],  ||1-t|| = %f',t,norm(1-t)),'fontsize',15)


VU13 = @(rho,m,E)rho.*(E.*rho-m.^2)./(E.*rho.*2.0-m.^2).^2.*-4.0;
VU23 = @(rho,m,E)m.*rho.^2.*1.0./(E.*rho.*2.0-m.^2).^2.*-4.0;
VU33 = @(rho,m,E)rho.^3.*1.0./(E.*rho.*2.0-m.^2).^2.*4.0;

clf
hold on
plot(rp,Vp*rho)
plot(rp,rhoavg - VU13(Vp*rho,Vp*m,Vp*E).*(Vp*rho-rhoavg),'.')
% plot(rp,mavg - VU23(Vp*rho,Vp*m,Vp*E).*(Vp*m-mavg),'.')
% plot(rp,Eavg - VU33(Vp*rho,Vp*m,Vp*E).*(Vp*E-Eavg),'.')

% plot(rp,VU13(Vp*rho,Vp*m,Vp*E),'.')
% plot(rp,VU23(Vp*rho,Vp*m,Vp*E),'.')
% plot(rp,VU33(Vp*rho,Vp*m,Vp*E),'.')

% VU13(rhoq,mq,Eq)
% VU23(rhoq,mq,Eq)
% VU33(rhoq,mq,Eq)
% VU13(Vf1*rho,Vf1*m,Vf1*E)
% VU23(Vf1*rho,Vf1*m,Vf1*E)
% VU33(Vf1*rho,Vf1*m,Vf1*E)
