N = 4;
r = JacobiGQ(0,0,N);
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

% % check min/max eigs
% A0 = zeros(3,3);
% for i = 1:10
%     rho = 10/2^i;
%     m = .25;
%     e = 10;
%     A0(1,1) = dUdV11(V1(rho,m,e),V2(rho,m,e),V3(rho,m,e));
%     A0(1,2) = dUdV12(V1(rho,m,e),V2(rho,m,e),V3(rho,m,e));
%     A0(1,3) = dUdV13(V1(rho,m,e),V2(rho,m,e),V3(rho,m,e));
%     A0(2,2) = dUdV22(V1(rho,m,e),V2(rho,m,e),V3(rho,m,e));
%     A0(2,3) = dUdV23(V1(rho,m,e),V2(rho,m,e),V3(rho,m,e));
%     A0(3,3) = dUdV33(V1(rho,m,e),V2(rho,m,e),V3(rho,m,e));
%     A0 = A0 + triu(A0)';
%     loglog(rho,min(eig(A0)),'o')
%     hold on
%     plot(rho,max(eig(A0)),'x')
% end
% return

k = pi;
rho = 2 + Pq*sin(k*rq);
m = Pq*sin(k*rq);
E = 2 + Pq*sin(k*rq);

% k = 0;
F0 = V*diag([1;zeros(N,1)])/V;

theta = 1;
rho = F0*rho + theta*(rho-F0*rho);
% m = F0*m + theta*(m-F0*m);
E = F0*E + theta*(E-F0*E);

rhoq = Vq*rho; 
mq = Vq*m;
Eq = Vq*E;
q1 = Pq*V1(rhoq,mq,Eq);
q2 = Pq*V2(rhoq,mq,Eq);
q3 = Pq*V3(rhoq,mq,Eq);


dUV11 = dUdV11(Vq*q1,Vq*q2,Vq*q3);
dUV12 = dUdV12(Vq*q1,Vq*q2,Vq*q3);
dUV13 = dUdV13(Vq*q1,Vq*q2,Vq*q3);
dUV22 = dUdV22(Vq*q1,Vq*q2,Vq*q3);
dUV23 = dUdV23(Vq*q1,Vq*q2,Vq*q3);
dUV33 = dUdV33(Vq*q1,Vq*q2,Vq*q3);
dq1 = V1(rhoq,mq,Eq)-Vq*q1;
dq2 = V2(rhoq,mq,Eq)-Vq*q2;
dq3 = V3(rhoq,mq,Eq)-Vq*q3;
dU1 = dUV11.*dq1 + dUV12.*dq2 + dUV13.*dq3;
dU2 = dUV12.*dq1 + dUV22.*dq2 + dUV23.*dq3;
dU3 = dUV13.*dq1 + dUV23.*dq2 + dUV33.*dq3;
dUv = [dU1';dU2';dU3']

rhovq = U1(Vq*q1,Vq*q2,Vq*q3);
mvq = U2(Vq*q1,Vq*q2,Vq*q3);
Evq = U3(Vq*q1,Vq*q2,Vq*q3);

dU1 = rhovq - rhoq;
dU2 = mvq - mq;
dU3 = Evq - Eq;
dU = -[dU1';dU2';dU3']
% dUv./dU

rhovp = U1(Vp*q1,Vp*q2,Vp*q3);
mvp = U2(Vp*q1,Vp*q2,Vp*q3);
Evp = U3(Vp*q1,Vp*q2,Vp*q3);

plot(rp,Vp*rho)
hold on
plot(rp,rhovp,'--');

figure
hold on
plot(rp,V3(Vp*rho,Vp*m,Vp*E))
plot(rp,Vp*q3,'--');
% plot(rp,ones(size(rp))*wq'*(Vq*q3)/2,'k-.') % mean of q3

rhoavg = wq'*(Vq*rho)/2;
mavg = wq'*(Vq*m)/2;
Eavg = wq'*(Vq*E)/2;

v3avg = wq'*V3(Vq*rho,Vq*m,Vq*E)/2;
plot(rp,ones(size(rp))*v3avg,'b--') % mean of q3
plot(rp,ones(size(rp))*V3(rhoavg,mavg,Eavg),'r-.') % mean of q3

return

plot(rp,Vp*m)
plot(rp,mvp,'--')

plot(rp,Vp*E)
plot(rp,Evp,'--')
xlim([-1.1 1.1])