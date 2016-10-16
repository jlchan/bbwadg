clear

N = 15;

r = JacobiGL(0,0,N);
[rq wq] = JacobiGQ(0,0,N+1);

re = linspace(-1,1,N+1)';
rp = linspace(-1,1,1000)';
VB = bern_basis_1D(N,r);
Vp = bern_basis_1D(N,rp);
VBe = bern_basis_1D(N,re);
Vq = bern_basis_1D(N,rq);

a = .25;
d = 1500;

f = @(x)  .5 + 1./(1 + exp(-d*(x-a)));%x > -1/3 & x < 1/3;%exp(-25*x.^2);
f = @(x) abs(x-a);
% f = @(x) sin(2*pi*x+.25);
u = VB\f(r);

plot(rp,f(rp),'-','linewidth',2)
hold on
plot(rp,Vp*u,'--','linewidth',2)

% plot(re,u,'o--');

TV = 0;
for i = 1:N
    TV = TV + abs(u(i,:) - u(i+1,:));
end
title(sprintf('TV = %f, max err = %g\n',TV,max(abs(f(rp)-Vp*u))))

W = get_BB_smoother(N);
% W = get_BB_P1smoother(N);
W = VBe;

% return
uB = W*u;
plot(rp,Vp*uB,'-.','linewidth',2);
% plot(re,uB,'o');

return
M = Vq'*diag(wq)*Vq;
uB = (M*u)./(Vq'*wq);
% uB = Vq'*(wq.*f(rq)) ./ (wq'*Vq)';
plot(rp,Vp*uB,'--');
plot(re,uB,'o');



hold off

sum(u)/(N+1)
sum(uB)/(N+1)
