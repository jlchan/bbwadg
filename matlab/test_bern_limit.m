N = 8;

r = JacobiGL(0,0,N);

V = bern_basis_1D(N,r);
re = linspace(-1,1,N+1)';
Ve = bern_basis_1D(N,re);

rp = linspace(-1,1,150)';
Vp = bern_basis_1D(N,rp);

f = @(x) 1./(1+exp(-100*(x-.1)));

err = max(abs(f(r) - Ve*f(re)))

a = err/N;
u = a*(V\f(r)) + f(r)*(1-a);
clf
plot(rp,f(rp),'-')
hold on
plot(re,Ve*u,'o')
plot(rp,Vp*u,'--')
drawnow
