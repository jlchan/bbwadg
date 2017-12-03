clear
N = 6;
r = JacobiGQ(0,0,N);
[rq wq] = JacobiGQ(0,0,N+4);
% rq = [-1;rq;1]; wq = [0;wq;0];
Nq = length(rq);

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

S = @(rho,m,E) -rho.*s(rho,m,E);

% conjugate entropy - Legendre transform
TV = @(V1,V2,V3) U1(V1,V2,V3).*V1 + U2(V1,V2,V3).*V2 + U3(V1,V2,V3).*V3 - S(U1(V1,V2,V3),U2(V1,V2,V3),U3(V1,V2,V3));
TU = @(U1,U2,U3) U1.*V1(U1,U2,U3) + U2.*V2(U1,U2,U3) + U3.*V3(U1,U2,U3) - S(U1,U2,U3);

% define vars
k = pi;
rhoex = @(x) 1.1 - exp(-100*(x-.25).^2);
mex = @(x) .1+sin(k*x);
Eex = @(x) 2 + sin(k*x);
rho = Pq*rhoex(rq);
m = Pq*mex(rq);
E = Pq*Eex(rq);

rhoavg = wq'*(Vq*rho)/2;
mavg = wq'*(Vq*m)/2;
Eavg = wq'*(Vq*E)/2;

rhop = Vp*rho;
mp = Vp*m;
Ep = Vp*E;

lim = @(u,theta) .5*wq'*u + theta*(u - .5*wq'*u);

tt = [0 1 0];
rhoq = lim(Vq*rho,tt(1));
mq = lim(Vq*m,tt(2));
Eq = lim(Vq*E,tt(3));

v1 = V1(rhoq,mq,Eq);
v2 = V2(rhoq,mq,Eq);
v3 = V3(rhoq,mq,Eq);
pv1 = Vq*Pq*v1;
pv2 = Vq*Pq*v2;
pv3 = Vq*Pq*v3;

% plot(rp,Vp*Pq*rhoq,'o-')
% hold on
% plot(rp,Vp*Pq*mq,'s-')
% plot(rp,Vp*Pq*Eq,'^-')
% legend('U1','U2','U3')
% return
% 
% plot(rp,V1(rhop,mp,Ep),'o-')
% hold on
% plot(rp,V2(rhop,mp,Ep),'s-')
% plot(rp,V3(rhop,mp,Ep),'^-')
% plot(rp,Vp*Pq*pv1,'.')
% plot(rp,Vp*Pq*pv2,'x')
% plot(rp,Vp*Pq*pv3,'*')
% legend('V1','V2','V3')
% return

Savg = -rhoavg.*s(rhoavg,mavg,Eavg);
Savg = 0;
Sq = -rhoq.*s(rhoq,mq,Eq)+Savg;

V1a = V1(rhoavg,mavg,Eavg);
V2a = V2(rhoavg,mavg,Eavg);
V3a = V3(rhoavg,mavg,Eavg);
U1v = @(rho,m,E) U1(Vq*Pq*V1(rho,m,E),Vq*Pq*V2(rho,m,E),Vq*Pq*V3(rho,m,E));
U2v = @(rho,m,E) U2(Vq*Pq*V1(rho,m,E),Vq*Pq*V2(rho,m,E),Vq*Pq*V3(rho,m,E));
U3v = @(rho,m,E) U3(Vq*Pq*V1(rho,m,E),Vq*Pq*V2(rho,m,E),Vq*Pq*V3(rho,m,E));

rhov = U1v(rhoq,mq,Eq);
mv = U2v(rhoq,mq,Eq);
Ev = U3v(rhoq,mq,Eq);
Svq = -rhov.*s(rhov,mv,Ev)+Savg;

% plot(rq,rhoex(rq),'o-')
% return
hold on
plot(rq,S(rhoq,mq,Eq),'o-')
plot(rq,S(rhov,mv,Ev),'x-')
% plot(rq,TU(rhoq,mq,Eq),'o-')
% plot(rq,TU(rhov,mv,Ev),'x-')
return

du1 = (rhoq-rhov);
du2 = (mq-mv);
du3 = (Eq-Ev);
dS = v1.*du1 + v2.*du2 + v3.*du3;
[~,iminS] = min(dS);
dS1 = v1.*du1;
dS2 = v2.*du2;
dS3 = v3.*du3;
dSi = [dS1(iminS) dS2(iminS) dS3(iminS)]; 

dv1 = (eye(size(Pq,2))-Vq*Pq)*V1(rhoq,mq,Eq);
dv2 = (eye(size(Pq,2))-Vq*Pq)*V2(rhoq,mq,Eq);
dv3 = (eye(size(Pq,2))-Vq*Pq)*V3(rhoq,mq,Eq);
% S(uv) - S(u) >
(Sq-Svq) > (pv1.*du1 + pv2.*du2 + pv3.*du3)
(Sq-Svq) < (v1.*du1 + v2.*du2 + v3.*du3) % want to guarantee Sq >= Svq

a = 2; % don't allow over (a)-times max increase in pointwise entropy 
Sratio = max(Svq)/max(Sq);
flag = Sratio > a; 

scale = (Sratio-a)/Sratio;
tvec = min(1, max(0,1 + (1-flag) + scale*dSi/max(abs(Svq))*flag));
%tvec = min(1, max(0,1 + (1-flag) + scale*dSi/max(abs(Svq))*flag));
% tvec(:) = .9;
% tvec = [1 1 1];
t1 = tvec(1); t2 = tvec(2); t3 = tvec(3);


rhov2 = U1v(lim(rhoq,t1),lim(mq,t1),lim(Eq,t3));
mv2   = U2v(lim(rhoq,t1),lim(mq,t1),lim(Eq,t3));
Ev2   = U3v(lim(rhoq,t1),lim(mq,t1),lim(Eq,t3));
Sv2q  = -rhov2.*s(rhov2,mv2,Ev2) + Savg;

% disp('average entropies: Savg, Sq, Svq, Sv2q')
% [Savg,wq'*Sq/2,wq'*Svq/2,wq'*Sv2q/2]
plot(rq,Sq,'o--')
hold on
plot(rq,Svq,'s--')
plot(rq,Sv2q,'x--')
% title(sprintf('Avg entropy: orig = %f, new = %f, limited = %f',wq'*Sq/2,wq'*Svq/2,wq'*Sv2q/2))
title(sprintf('Max entropy: orig = %f, new = %f, lim = %f',max(Sq),max(Svq),max(Sv2q)))
% plot(rp,Savg*ones(size(rp)))
return

figure
hold on
plot(rp,Vp*rho,'-')
plot(rq,rhov,'o');
plot(rp,U1(Vp*Pq*V1(rhoq,mq,Eq),Vp*Pq*V2(rhoq,mq,Eq),Vp*Pq*V3(rhoq,mq,Eq)),'--');

rhoq = Vq*(rhoavg + t1*(rho-rhoavg));
mq = Vq*(mavg + t2*(m-mavg));
Eq = Vq*(Eavg + t3*(E-Eavg));
Sq2 = -rhoq.*s(rhoq,mq,Eq)+Savg;
plot(rp,U1(Vp*Pq*V1(rhoq,mq,Eq),Vp*Pq*V2(rhoq,mq,Eq),Vp*Pq*V3(rhoq,mq,Eq)),'-.');

disp('max ptwise entropies: Sq, Svq, Sq2, Sv2q')
[max(Sq),max(Svq),max(Sq2),max(Sv2q)]

