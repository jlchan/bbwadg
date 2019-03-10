clear

N = 4;
% Ntarget = N+1;
Nq = (N+1);
rqref = linspace(-1,1,10000)';
wqref = [.5; ones(length(rqref)-2,1); .5]; wqref = wqref/sum(wqref)*2;
% [rq wq] = JacobiGQ(0,0,Nq);
Vqref = Vandermonde1D(Nq,rqref);

VDEIM(:,1) = Vqref(:,1);
[~,id] = max(abs(Vqref(:,1)));
p = id;
for j = 2:Nq+1
    r = Vqref(:,j)-VDEIM*(VDEIM(p,:)\Vqref(p,j));
    [~,id] = max(abs(r));
    p(j) = id;
    VDEIM = [VDEIM r];    
end
rq = rqref(p);

% use as quad pts
Vq = Vandermonde1D(N,rq);
Vf = Vandermonde1D(N,[-1;1]);
D = Vandermonde1D(N,JacobiGL(0,0,N))\GradVandermonde1D(N,JacobiGL(0,0,N));
wf = [1;1];

% min || VDEIM(p,:)'*w - VDEIM'*wq ||
% constraint C*w = d
A = Vqref(p,:)*Vqref(p,:)';
b = Vqref(p,:)*Vqref'*wqref;
C = D'*Vq';
d = Vf'*diag([-1;1])*wf;

S = [A C'
 C zeros(size(C,1))];
W = pinv(S)*[b;d];
wq = W(1:Nq+1);

M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq));
E = Vf*Pq;
B = diag([-1;1]);
Q = diag(wq)*Vq*D*Pq;

e = ones(size(Q,1),1);
norm(e'*Q - e'*E'*B*E)

plot(rq,wq,'o')
hold on
plot([-1,1],[0 0],'--')





