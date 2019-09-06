clear

N = 50;
% Ntarget = N+1;
Nq = (N+1);
rqref = linspace(-1,1,10000)';
wqref = [.5; ones(length(rqref)-2,1); .5]; wqref = wqref/sum(wqref)*2;
% [rq wq] = JacobiGQ(0,0,Nq);
Vqref = Vandermonde1D(Nq,rqref);

id = get_empirical_cubature(Vqref,wqref,1e-7,inf);
rq = rqref(id);
[rq,p] = sort(rq);
id = id(p);
wq0 = Vqref(id,:)' \ (Vqref'*wqref);

% use as quad pts
Vq = Vandermonde1D(N,rq);
Vf = Vandermonde1D(N,[-1;1]);
D = Vandermonde1D(N,JacobiGL(0,0,N))\GradVandermonde1D(N,JacobiGL(0,0,N));
wf = [1;1];

% min || V(id,:)'*w - V'*wq ||
% constraint C*w = d 
%        ->  Q'*1 = E'*B*1
A = Vqref(id,:)*Vqref(id,:)';
b = Vqref(id,:)*Vqref'*wqref;
C = D'*Vq';
d = Vf'*diag([-1;1])*wf;

S = [A C'
 C zeros(size(C,1))];
W = pinv(S)*[b;d];
wq = W(1:Nq+1);
% wq = wq0;

M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq));
E = Vf*Pq;
B = diag([-1;1]);
Q = diag(wq)*Vq*D*Pq;

e = ones(size(Q,1),1);
norm(e'*Q - e'*E'*B*E)

plot(rq,wq,'o')
hold on
plot(rq,wq0,'x')
title(sprintf('diff = %g\n',norm(wq-wq0)))
plot([-1,1],[0 0],'--')





