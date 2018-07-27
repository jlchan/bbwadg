clear
N = 7;
r = JacobiGL(0,0,N);
[rq wq] = JacobiGQ(0,0,N+1);
% rq = [-1;rq;1]; wq = [0;wq;0];
Nq = length(rq);

V = Vandermonde1D(N,r);
Vq = Vandermonde1D(N,rq)/V;
Vf = Vandermonde1D(N,[-1;1])/V;
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


VU = @(U) [V1(U(:,1),U(:,2),U(:,3)),V2(U(:,1),U(:,2),U(:,3)),V3(U(:,1),U(:,2),U(:,3))];
UV = @(V) [U1(V(:,1),V(:,2),V(:,3)),U2(V(:,1),V(:,2),V(:,3)),U3(V(:,1),V(:,2),V(:,3))];

U1v = @(V,rho,m,E) U1(V*Pq*V1(rho,m,E),V*Pq*V2(rho,m,E),V*Pq*V3(rho,m,E));
U2v = @(V,rho,m,E) U2(V*Pq*V1(rho,m,E),V*Pq*V2(rho,m,E),V*Pq*V3(rho,m,E));
U3v = @(V,rho,m,E) U3(V*Pq*V1(rho,m,E),V*Pq*V2(rho,m,E),V*Pq*V3(rho,m,E));

% define vars
a = 1;
b = 0;
k = 1;
%rho = 1 - b*Pq*exp(-10*(rq).^2);
% rho = 1 - b*Pq*sin(a+k*(rq));
rho = 1.000 - b*Pq*rq;
m = 0 + b*Pq*sin(a+k*rq);
E = 1 + b*Pq*exp(-10*(rq).^2);
% m = 1 + 0*r;
% E = 1 + 0*r;

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

VU13 = @(rho,m,E) rho.*(E.*rho-m.^2).*1.0./(E.*rho.*2.0-m.^2).^2.*-4.0;
VU23 = @(rho,m,E) -m.*1.0./(E-(m.^2.*(1.0./2.0))./rho).^2;
VU33 = @(rho,m,E) rho.*1.0./(E-(m.^2.*(1.0./2.0))./rho).^2;
dVU13 = @(U) VU13(U(:,1),U(:,2),U(:,3));
dVU23 = @(U) VU23(U(:,1),U(:,2),U(:,3));
dVU33 = @(U) VU33(U(:,1),U(:,2),U(:,3));


%%
% plot 
lim = @(u,theta) .5*wq'*u + theta.*(u - .5*wq'*u);
% s = @(rho,m,E) log((gamma-1).*(E - .5*m.^2./rho)./(rho.^gamma));
Sfun = @(rho,m,E) -rho.*s(rho,m,E);

SU = @(U) Sfun(U(:,1),U(:,2),U(:,3));

VV = @(U) [V1(U(:,1),U(:,2),U(:,3)),V2(U(:,1),U(:,2),U(:,3)),V3(U(:,1),U(:,2),U(:,3))];
UU = @(V) [U1(V(:,1),V(:,2),V(:,3)),U2(V(:,1),V(:,2),V(:,3)),U3(V(:,1),V(:,2),V(:,3))];
U = [rho m E];
plot(rp,Vp*U,'b--')
hold on
plot(rp,SU(Vp*U),'bo')
hold on
% plot(rp,SU(UU(Vp*Pq*VV(Vq*U))),'rx')

% plot(rp,SU(UU(Vp*Pq*VV(lim(Vq*U,.95)))),'k^')


return

%%
clf
hold on
plot(rp,V3(Vp*rho,Vp*m,Vp*E))
% [VU13(rhoavg,mavg,Eavg) VU23(rhoavg,mavg,Eavg) VU33(rhoavg,mavg,Eavg)]
% dV3 = VU13(Vp*rho,Vp*m,Vp*E).*(Vp*rho-rhoavg) + ...
%     VU23(Vp*rho,Vp*m,Vp*E).*(Vp*m-mavg) + ...
%     VU33(Vp*rho,Vp*m,Vp*E).*(Vp*E-Eavg);
% plot(rp,V3(rhoavg,mavg,Eavg) + dV3,'--')

plot(rp,Vp*Pq*V3(Vq*rho,Vq*m,Vq*E),'.-')

t = [1 1 1];
rholim = rhoavg + t(1)*(rho-rhoavg);
mlim = mavg + t(2)*(m-mavg);
Elim = Eavg + t(3)*(E-Eavg);

% k = 5;
% F = V*diag([ones((N-k),1);zeros(k+1,1)])*inv(V);
% rholim = F*rho;
% mlim = F*m;
% Elim = F*E;
plot(rp,Vp*Pq*V3(Vq*rholim,Vq*mlim,Vq*Elim),'x-')
legend('exact','proj')

axis([-1,1 -1.5 -.5])

%%
figure
clf
plot(rp,Vp*rho)
hold on
plot(rp,U1(Vp*Pq*V1(Vq*rho,Vq*m,Vq*E),Vp*Pq*V2(Vq*rho,Vq*m,Vq*E),Vp*Pq*V3(Vq*rho,Vq*m,Vq*E)),'--')

%% ================================ gradient descent limiting
e = [0 0 1]';
theta = .5*[1 1 1];
U = [rho m E];
Uq = Vq*U;
Uavg = .5*wq'*Uq;

for i = 1:1000
    Utheta = Uavg + (Uq-Uavg)*diag(theta);
    rt = Vf*Pq*(VU(Utheta)*e) - VU(Vf*U)*e;
    drt1 = Vf*Pq*(dVU13(Utheta));
    
    dt1 = rt'*drt1;
    dtheta = rt'*Vf*Pq*[(dVU13(Utheta)) (dVU23(Utheta)) (dVU33(Utheta))];
    alpha = .001;
    theta = theta - alpha*dtheta;
end
