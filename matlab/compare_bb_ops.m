clear

N = 10;

r = JacobiGL(0,0,N);
[rq wq] = JacobiGQ(0,0,N+1);

re = linspace(-1,1,N+1)';
rp = linspace(-1,1,200)';
VB = bern_basis_1D(N,r);
Vp = bern_basis_1D(N,rp);
VBe = bern_basis_1D(N,re);
Vq = bern_basis_1D(N,rq);

a = -.1;
d = 150;

% f = @(x)  .5*x + 1./(1 + exp(-d*(x-a)));%x > -1/3 & x < 1/3;%exp(-25*x.^2);
f = @(x) x;
% f = @(x) sin(pi*x+.25);
u = VB\f(r);

plot(rp,f(rp),'-')
hold on
plot(rp,Vp*u,'-')

% uB = VBe*u;
% plot(rp,Vp*uB,'.-');
% plot(re,uB,'o');
% return
% W = get_BB_smoother(N);
W = get_BB_P1smoother(N);

% return
uB = W*u;
plot(rp,Vp*uB,'-.');
plot(re,uB,'o');

return
M = Vq'*diag(wq)*Vq;
uB = (M*u)./(Vq'*wq);
% uB = Vq'*(wq.*f(rq)) ./ (wq'*Vq)';
plot(rp,Vp*uB,'--');
plot(re,uB,'o');



hold off

sum(u)/(N+1)
sum(uB)/(N+1)
