clear
N = 7;
r = JacobiGQ(0,0,N);
[rq wq] = JacobiGQ(0,0,N+1);
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
a = 1;
k = pi;
rho = 1.475 + Pq*sin(a+k*rq);
m = .1+Pq*sin(a+k*rq);
E = 2 + Pq*sin(a+k*rq);

rhoavg = wq'*(Vq*rho)/2;
mavg = wq'*(Vq*m)/2;
Eavg = wq'*(Vq*E)/2;

rhop = Vp*rho;
mp = Vp*m;
Ep = Vp*E;

lim = @(u,theta) .5*wq'*u + theta*(u - .5*wq'*u);

rhoq = Vq*rho; 
mq = Vq*m;
Eq = Vq*E;

% rhoq = lim(rhoq,.5);
% mq = lim(mq,1);
% Eq = lim(Eq,.5);

Savg = 0;
Sq = -rhoq.*s(rhoq,mq,Eq)+Savg;

V1a = V1(rhoavg,mavg,Eavg);
V2a = V2(rhoavg,mavg,Eavg);
V3a = V3(rhoavg,mavg,Eavg);
v1 = V1(rhoq,mq,Eq);
v2 = V2(rhoq,mq,Eq);
v3 = V3(rhoq,mq,Eq);
pv1 = Vq*Pq*v1;
pv2 = Vq*Pq*v2;
pv3 = Vq*Pq*v3;

U1v = @(V,rho,m,E) U1(V*Pq*V1(rho,m,E),V*Pq*V2(rho,m,E),V*Pq*V3(rho,m,E));
U2v = @(V,rho,m,E) U2(V*Pq*V1(rho,m,E),V*Pq*V2(rho,m,E),V*Pq*V3(rho,m,E));
U3v = @(V,rho,m,E) U3(V*Pq*V1(rho,m,E),V*Pq*V2(rho,m,E),V*Pq*V3(rho,m,E));

rhov = U1v(Vq,rhoq,mq,Eq);
mv = U2v(Vq,rhoq,mq,Eq);
Ev = U3v(Vq,rhoq,mq,Eq);
Svq = -rhov.*s(rhov,mv,Ev)+Savg;

% Jacobian on time derivative
dVdU11 = @(rho,m,E) (1.0./(E.*rho.*2.0-m.^2).^2.*(gamma.*m.^4+m.^4+E.^2.*gamma.*rho.^2.*4.0-E.*gamma.*m.^2.*rho.*4.0))./rho;
dVdU12 = @(rho,m,E) m.^3.*1.0./(E.*rho.*2.0-m.^2).^2.*-2.0;
dVdU13 = @(rho,m,E) rho.*(E.*rho-m.^2).*1.0./(E.*rho.*2.0-m.^2).^2.*-4.0;
dVdU22 = @(rho,m,E) rho.*(E.*rho.*2.0+m.^2).*1.0./(E.*rho.*2.0-m.^2).^2.*2.0;
dVdU23 = @(rho,m,E) m.*rho.^2.*1.0./(E.*rho.*2.0-m.^2).^2.*-4.0;
dVdU33 = @(rho,m,E) rho.^3.*1.0./(E.*rho.*2.0-m.^2).^2.*4.0;

HU = @(rho,m,E) V1(rho,m,E).*rho + V2(rho,m,E).*m + V3(rho,m,E).*E;

if 1
    plot(rp,Vp*rho,'-')
    hold on
%     plot(rp,Vp*m,'-')
%     plot(rp,Vp*E,'-')
    %     plot(rp,U1v(Vp,rhoq,mq,Eq),'--')
    %     plot(rp,U2v(Vp,rhoq,mq,Eq),'--')
    %     plot(rp,U3v(Vp,rhoq,mq,Eq),'--')
    plot(rp,Vp*Pq*U1v(Vq,rhoq,mq,Eq),'--')
%     plot(rp,Vp*Pq*U2v(Vq,rhoq,mq,Eq),'--')
%     plot(rp,Vp*Pq*U3v(Vq,rhoq,mq,Eq),'--')
    
    % plot(rp,V1(Vp*rho,Vp*m,Vp*E),'-')
    % hold on
    % plot(rp,V2(Vp*rho,Vp*m,Vp*E),'-')
    % plot(rp,V3(Vp*rho,Vp*m,Vp*E),'-')
    % plot(rp,Vp*Pq*V1(Vq*rho,Vq*m,Vq*E),'--')
    % plot(rp,Vp*Pq*V2(Vq*rho,Vq*m,Vq*E),'--')
    % plot(rp,Vp*Pq*V3(Vq*rho,Vq*m,Vq*E),'--')
    
    plot(rq,HU(rhoq,mq,Eq),'o-')
    hold on
    plot(rq,HU(rhov,mv,Ev),'x-')
    return
end

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

a = 8; % don't allow over (a)-times max increase in pointwise entropy 
Sratio = max(Svq)/max(Sq);
flag = Sratio > a; 

scale = (Sratio-a)/Sratio;
tvec = min(1, max(0,1 + (1-flag) + scale*dSi/max(abs(Svq))*flag));
t1 = tvec(1); t2 = tvec(2); t3 = tvec(3);

rhov2 = U1v(Vq,lim(rhoq,t1),lim(mq,t1),lim(Eq,t3));
mv2   = U2v(Vq,lim(rhoq,t1),lim(mq,t1),lim(Eq,t3));
Ev2   = U3v(Vq,lim(rhoq,t1),lim(mq,t1),lim(Eq,t3));
Sv2q  = -rhov2.*s(rhov2,mv2,Ev2) + Savg;

plot(rq,Sq,'o--')
hold on
plot(rq,Svq,'s--')
plot(rq,Sv2q,'x--')
title(sprintf('Max entropy: orig = %f, new = %f, lim = %f',max(Sq),max(Svq),max(Sv2q)))



