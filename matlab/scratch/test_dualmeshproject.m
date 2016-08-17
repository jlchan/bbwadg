N = 3;

r = JacobiGL(0,0,N);
V = bern_basis_1D(N,r);

u = randn(N+1,1);
% u = V\sin(r+1);

re = linspace(-1,1,N+1)';
rp = linspace(-1,1,250)';
Vp = bern_basis_1D(N,rp);

hold on
% plot(re,u,'o-'); plot(rp,Vp*u,'--')

u1 = randn(N+1,1); u2 = randn(N+1,1); hold on
% u1 = Vsub*u; u2 =rot90(Vsub,2)*u;
Vsub = V\bern_basis_1D(N,(r-1)/2);

rp1 = (rp-1)/2; rp2 = (1+rp)/2;
plot((re-1)/2,u1,'s'); plot((re+1)/2,u2,'s')
plot(rp1,Vp*u1); plot(rp2,Vp*u2)

% project onto single element
[rq1 w] = JacobiGQ(0,0,N);
rq = [(rq1-1)/2; (rq1+1)/2]; w = [w; w]/2; % 2-elem join

Vq = bern_basis_1D(N,rq);
M = Vq'*diag(w)*Vq;

Vq1 = bern_basis_1D(N,rq1);
join = M\(Vq'*diag(w)*blkdiag(Vq1,Vq1));

ujoin = join*[u1;u2];
norm(u-ujoin)
plot(rp,Vp*ujoin)