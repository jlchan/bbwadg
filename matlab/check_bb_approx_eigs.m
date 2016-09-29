clear
N = 7;

r = JacobiGL(0,0,N);
re = linspace(-1,1,N+1)';
rp = linspace(-1,1,200)';
[VB Vr] = bern_basis_1D(N,r);
Dr = VB\Vr;
VBe = bern_basis_1D(N,re);
Vp = bern_basis_1D(N,rp);
[rq wq] = JacobiGQ(0,0,N);
Vq = bern_basis_1D(N,rq);
M = Vq'*diag(wq)*Vq;

% VBe = get_BB_smoother(N);
[W D] = eig(VBe);
lam = diag(D);

[lam p] = sort(lam);
W = W(:,p);

% Vandermonde1D(N,r)\(VB*W)
A = W(:,1:N-1);
a = A'*M*A;
a(abs(a)<1e-8) = 0;

% plot(lam,'o-');  return

wp = Vp*(W);
wp = wp*diag(1./max(abs(wp),[],1));
plot(rp,wp(:,1:N+1),'-')

hold on

% Vp = Vandermonde1D(N,rp);
% Vp = Vp * diag (1./Vp(1,:));
% plot(rp,Vp(:,3:N),'.')

