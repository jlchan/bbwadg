N = 3;
r = JacobiGL(0,0,N);
[rq wq] = JacobiGQ(0,0,N+2);
% rq = [-1;rq;1]; wq = [0;wq;0];

V = Vandermonde1D(N,r);
Vq = Vandermonde1D(N,rq)/V;
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));

% Vq = eye(N+1);
% Pq = Vq;

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
plot(rp,Vp*u)
plot(rp,uv,'--')

plot(rp,Vp*v)
plot(rp,vfun(Vp*u),'--')

% run one step of Newton for max, one for min - can do simultaneously?
[vmax id] = max(vfun(Vq*u));
% [vmax id] = max(vfun(u));

uavg = wq'*(Vq*u)/2;

% how to choose initial theta?
theta = .5;
for i = 1:1
    utheta = uavg + theta*(u-uavg);
    dfdt = Pq*(dvdu(Vq*utheta).*(Vq*(u-uavg)));
%     dfdt = (dvdu(utheta).*((u-uavg)));
    fdf = Vq*(Pq*vfun(Vq*utheta) - vmax)./(Vq*dfdt);
%     fdf = (Pq*vfun(Vq*utheta) - vmax)./(dfdt);
    theta = theta - fdf(id);
end
theta
utheta = uavg + theta*(u-uavg);

plot(rp,Vp*Pq*vfun(Vq*utheta),'.')
plot(rp,ufun(Vp*Pq*vfun(Vq*utheta)),'.')



