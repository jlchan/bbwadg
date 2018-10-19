N = 13;
[r s] = Nodes2D(N); [r s] = xytors(r,s);
[rq sq wq] = Cubature2D(2*N);
V = Vandermonde2D(N,r,s);
Vq = Vandermonde2D(N,rq,sq)/V;
[Vr Vs] = GradVandermonde2D(N,r,s);
Dr = Vr/V; Ds = Vs/V;

Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));

Drq = Vq*Dr*Pq;

f = @(x,y) exp(x+y);
dfdx = @(x,y) exp(x+y);
norm(Vq*Pq*f(rq,sq) - f(rq,sq),'fro')
norm(Drq*f(rq,sq) - dfdx(rq,sq),'fro')

%%

[a1D,wa] = JacobiGQ(0,0, N);
[b1D,wb] = JacobiGQ(1,0, N);
a = ones(length(a1D),1)*a1D';
b = b1D*ones(1,length(a1D));

rq = 0.5*(1+a).*(1-b)-1;
sq = b;
wq = 0.5*wb*(wa');

rq = rq(:); sq = sq(:); wq = wq(:);

aa = 2*(1+rq)./(1-sq)-1;
bb = sq;
% plot(aa,bb,'o')

% dudr = dr/da*du/da + dr/db * du/db
drda = (1-bb)/2; 
drdb = -(1+aa)/2;
Da = GradVandermonde1D(N,a1D)/Vandermonde1D(N,a1D);
Db = GradVandermonde1D(N,b1D)/Vandermonde1D(N,b1D);

Da = kron(Da,eye(N+1));
Db = kron(eye(N+1),Db);

Dr = diag(drda)*Da + diag(drdb)*Db;
Ds = Db; 

Vq = Vandermonde2D(N,rq,sq)/V;
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));

f = @(x,y) exp(x+y);
dfdx = @(x,y) exp(x+y);

norm(Vq*Pq*f(rq,sq) - f(rq,sq),'fro')
norm(Vq*Pq*(Ds*f(rq,sq)) - f(rq,sq),'fro')

vv = dfdx(rq,sq);
color_line3(aa,bb,vv,vv,'o')

vv = Vq*Pq*(Ds*f(rq,sq));
color_line3(aa,bb,vv,vv,'.')
