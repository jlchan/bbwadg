clear
N = 4;
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
k = pi;
rho = 1 - .4*Pq*exp(-10*rq.^2);
m = 1 + Pq*sin(a+k*rq);
E = 2 + Pq*sin(a+k*rq);
% rho = 1+0*r;
% m = 1+0*r; 
% E = 2+0*r;

rhoavg = wq'*(Vq*rho)/2;
mavg = wq'*(Vq*m)/2;
Eavg = wq'*(Vq*E)/2;

rhoq = Vq*rho;
mq = Vq*m;
Eq = Vq*E;

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

% minimize .5*(S(uv)-S(u))^2 + (t-1).^2
t = [1 1 1];
Vf1 = Vandermonde1D(N,1)/V;

S = @(rho,m,E) -rho.*s(rho,m,E);
% conjugate entropies 
TV = @(V1,V2,V3) U1(V1,V2,V3).*V1 + U2(V1,V2,V3).*V2 + U3(V1,V2,V3).*V3 - S(U1(V1,V2,V3),U2(V1,V2,V3),U3(V1,V2,V3));
TU = @(U1,U2,U3) U1.*V1(U1,U2,U3) + U2.*V2(U1,U2,U3) + U3.*V3(U1,U2,U3) - S(U1,U2,U3);
TTV = @(V1,V2,V3) U1(V1,V2,V3).*V1 + U2(V1,V2,V3).*V2 + U3(V1,V2,V3).*V3;
TTU = @(U1,U2,U3) U1.*V1(U1,U2,U3) + U2.*V2(U1,U2,U3) + U3.*V3(U1,U2,U3);

S0 = max(S(Vq*rho,Vq*m,Vq*E));
J = @(t) (S(U1v(Vf1,lim(rhoq,t(1)),lim(mq,t(2)),lim(Eq,t(3))),...
    U2v(Vf1,lim(rhoq,t(1)),lim(mq,t(2)),lim(Eq,t(3))),...
    U3v(Vf1,lim(rhoq,t(1)),lim(mq,t(2)),lim(Eq,t(3)))) - S0)^2;

% T0 = max(TU(Vq*rho,Vq*m,Vq*E));
% J = @(t) TV(Vf1*Pq*V1(lim(rhoq,t(1)),lim(mq,t(2)),lim(Eq,t(3))),...
%     Vf1*Pq*V2(lim(rhoq,t(1)),lim(mq,t(2)),lim(Eq,t(3))),...
%     Vf1*Pq*V3(lim(rhoq,t(1)),lim(mq,t(2)),lim(Eq,t(3)))- T0)^2;

% T0 = max(TTU(Vq*rho,Vq*m,Vq*E));
% J = @(t) TTV(Vf1*Pq*V1(lim(rhoq,t(1)),lim(mq,t(2)),lim(Eq,t(3))),...
%     Vf1*Pq*V2(lim(rhoq,t(1)),lim(mq,t(2)),lim(Eq,t(3))),...
%     Vf1*Pq*V3(lim(rhoq,t(1)),lim(mq,t(2)),lim(Eq,t(3)))- T0)^2;

% [t,fval,flag,out] = fminsearch(@(t) J(t) + .1*norm(t-1)^2,t);
t = .9*t;
rhoq = lim(rhoq,t(1));
mq = lim(mq,t(2));
Eq = lim(Eq,t(3));

du1 = rhoq-rhovq;
du2 = mq-mvq;
du3 = Eq-Evq;
dv1 = v1-pv1;
dv2 = v2-pv2;
dv3 = v3-pv3;

dT = rhovq.*dv1 + mq.*dv2 + Evq.*dv3;

% dU = [max(abs(du1)),max(abs(du2)),max(abs(du3))]; dU = dU/max([max(abs(pv1)),max(abs(pv2)),max(abs(pv3))])
% dV = [max(abs(dv1)),max(abs(dv2)),max(abs(dv3))]; dV = dV/max([max(abs(rhovq)),max(abs(mvq)),max(abs(Evq))])

hold on
% plot(rp,Vp*rho,'-','linewidth',2)
% plot(rp,Vp*m,'-')
% plot(rp,Vp*E,'-')
% plot(rp,U1v(Vp,rhoq,mq,Eq),'--','linewidth',2)
% plot(rp,U2v(Vp,rhoq,mq,Eq),'--')
% plot(rp,U3v(Vp,rhoq,mq,Eq),'--')
% plot(rp,Vp*Pq*U1v(Vq,rhoq,mq,Eq),'--')
% plot(rp,Vp*Pq*U2v(Vq,rhoq,mq,Eq),'--')
% plot(rp,Vp*Pq*U3v(Vq,rhoq,mq,Eq),'--')
% plot(rp,S(U1v(Vp,Vq*rho,Vq*m,Vq*E),U2v(Vp,Vq*rho,Vq*m,Vq*E),U3v(Vp,Vq*rho,Vq*m,Vq*E)),'.-','linewidth',2)
plot(rp,S(Vp*rho,Vp*m,Vp*E),'-','linewidth',2)
plot(rq,S(Vq*rho,Vq*m,Vq*E),'o','linewidth',2)
plot(rp,ones(size(rp))*S0,'linewidth',2)
plot(rp,S(U1v(Vp,rhoq,mq,Eq),U2v(Vp,rhoq,mq,Eq),U3v(Vp,rhoq,mq,Eq)),'--','linewidth',2)
title(sprintf('t = [%1.2f %1.2f %1.2f],  ||1-t|| = %f',t,norm(1-t)),'fontsize',15)
