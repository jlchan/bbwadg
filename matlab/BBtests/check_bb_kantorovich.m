clear

Globals1D
N = 6;

K1D = N+1; % macro grid
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
% VX = JacobiGL(0,0,N+1)';

% Initialize solver and construct grid and metric
StartUp1D;

[rq wq] = JacobiGQ(0,0,N+1);
Vq = Vandermonde1D(N,rq)/V;
xq = Vq*x;

re = linspace(-1,1,N+1)';
rp = linspace(-1,1,200)';
xp = Vandermonde1D(N,rp)/V * x;
VB = bern_basis_1D(N,r);
Vp = bern_basis_1D(N,rp);
VBe = bern_basis_1D(N,re);

V_sub = bern_basis_1D(N,x(:)); % interpolate to micro grid

% f = @(x) sin(pi*x);
a = .1;
% avec = -.5:.1:1;
% for i = 1:length(avec)
%     a = avec(i);
d = 50;

f = @(x) 1 + 1./(1 + exp(-d*(x-a)));%x > -1/3 & x < 1/3;%exp(-25*x.^2);
% f = @(x) sin(pi*(x-a));
u = VB\f(r);

plot(rp,f(rp),'-')
hold on
plot(rp,Vp*u,'-')

%     plot(xq,0*xq,'x');

wJ = diag(wq)*(Vq*J);

uB = (N+1)/2 * sum(wJ.*(Vq*reshape(V_sub*u,N+1,K)),1);
uB = uB(:);
%     uB = f(re);

plot(rp,Vp*uB,'--');
plot(re,uB,'o');

hold off
% end
% check to make sure conservation is maintained
% sum(u,1)/(N+1)
% sum(uB,1)/(N+1)

utmp = zeros(N+1,1);
for i = 1:N+1
    utmp(i) = 1;
    uB = (N+1)/2 * sum(wJ.*(Vq*reshape(V_sub*utmp,N+1,K)),1);
    W(:,i) = uB(:);
    utmp(i) = 0;
end
% plot(rp,Vp*We,'-'); return

clf
plot(rp,f(rp),'-')
hold on
plot(rp,Vp*u,'-')
plot(rp,Vp*VBe*u,'.-')
plot(rp,Vp*W*u,'--');
plot(re,W*u,'o');

sum(W*u)/(N+1)
sum(u)/(N+1)
