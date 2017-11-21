N = 3;
r = JacobiGL(0,0,N);
[rq wq] = JacobiGQ(0,0,N+1);
rq = [-1;rq;1]; wq = [0;wq;0];

V = Vandermonde1D(N,r);
Vq = Vandermonde1D(N,rq)/V;
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));

rp = linspace(-1,1,250)';
Vp = Vandermonde1D(N,rp)/V;

f = @(x) 1.01 + sin(pi*x);
u = Pq*f(rq);

vfun = @(u) log(u);
dvdu = @(u) 1./(u);
ufun = @(v) exp(v);
dudv = @(v) exp(v);

v = Pq*vfun(Vq*u);
uv = ufun(Vp*v);

hold on
plot(rp,Vp*u0)
plot(rp,uv,'--')

plot(rp,Vp*v)
plot(rp,vfun(Vp*u0),'--')

[val id] = max(uv);
[vmax id] = max(vfun(Vp*u));

uavg = wq'*(Vq*u)/2;

% how to choose initial theta?
theta = .5;
utheta = uavg + theta*(u-uavg);
dfdt = Pq*(dvdu(Vq*utheta).*(Vq*(u-uavg)));

fdf = (Vp*Pq*vfun(Vq*utheta) - vmax)./(Vp*dfdt);
theta = theta - fdf(id);
theta
utheta = uavg + theta*(u-uavg);

plot(rp,Vp*Pq*vfun(Vq*utheta),'.')
plot(rp,ufun(Vp*Pq*vfun(Vq*utheta)),'.')
