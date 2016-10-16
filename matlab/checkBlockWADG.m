N = 7;
r = JacobiGL(0,0,N);
[rq wq] = JacobiGQ(0,0,N+1);

V = Vandermonde1D(N,r);
Vq = Vandermonde1D(N,rq)/V;
invM = V*V';
M = inv(invM);

w = @(x) 1 + .1*sin(pi*x);

Mw = Vq'*diag(wq.*w(rq))*Vq;
Mwadg = M*(Vq'*diag(wq./w(rq))*Vq)\M;
invMwadg = invM*(Vq'*diag(wq./w(rq))*Vq)*invM;
b = Vq'*(wq.*exp(rq));
u1 = Mw\b;
u2 = invMwadg*b;

norm(u1-u2)

%% coupled mass matrix

a = @(x) 2 + .1*sin(pi*x);
b = @(x) 1 - .1*sin(pi*x);
c = @(x) 2 + .1*sin(pi*x);

det = @(x) a(x).*c(x)-b(x).^2;
ia = @(x) c(x)./det(x);
ib = @(x) -b(x)./det(x);
ic = @(x) a(x)./det(x); 

Ma = Vq'*diag(wq.*a(rq))*Vq;
Mb = Vq'*diag(wq.*b(rq))*Vq;
Mc = Vq'*diag(wq.*c(rq))*Vq;
MW = [Ma Mb; Mb Mc];

Mia = Vq'*diag(wq.*ia(rq))*Vq;
Mib = Vq'*diag(wq.*ib(rq))*Vq;
Mic = Vq'*diag(wq.*ic(rq))*Vq;

MWadg = blkdiag(invM,invM) * [Mia Mib; Mib Mic] * blkdiag(invM,invM);
b = Vq'*(wq.*exp(rq));
b = [b;2*b];

u1 = MW\b;
u2 = MWadg*b;
norm(u1-u2)

