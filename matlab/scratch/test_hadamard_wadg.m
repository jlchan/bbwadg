% clc
N = 3;

r = JacobiGL(0,0,N);
% r = JacobiGQ(0,0,N);

% [rq wq] = JacobiGL(0,0,N); rq = .9999999999*rq;
[rq wq] = JacobiGQ(0,0,N+1);
Nq = length(rq);

V = Vandermonde1D(N,r);
Dr = GradVandermonde1D(N,r)/V;
Vq = Vandermonde1D(N,rq)/V;

M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq));

w = 2+sin(1+pi*rq);
% w = Vq*Pq*w;
uq = Vq*r;
u = Pq*uq;

Minvw = Vq'*diag(wq./w)*Vq;
Tinvw = Vq*(Minvw\(Vq'*diag(wq)));

% wq'*(w.*uq)
% wq'*(Tinvw*uq)

% wu = repmat(uq',N+1,1).*repmat((1./w)',N+1,1);
% sum(Pq.*wu,2)

Pq*diag(1./w)*uq
diag(Pq*diag(1./w)*Vq*repmat(u,1,N+1))
(Pq.*(Vq*repmat(u,1,N+1))')*(1./w)
