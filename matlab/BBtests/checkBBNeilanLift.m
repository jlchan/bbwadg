clear
N = 4;


[rfq wfq] = JacobiGQ(0,0,N+1);
sfq = -ones(size(rfq));
[rq sq wq] = Cubature2D(2*N+1);
Vq = bern_basis_tri(N,rq,sq);
Vfq = bern_basis_tri(N,rfq,sfq);
M = Vq'*diag(wq)*Vq;
Mf = Vfq'*diag(wfq)*Vfq;

LIFT = M\Mf(:,1:N+1);
LIFT(abs(LIFT)<1e-8) = 0;

% increased order lift
Vfq2 = bern_basis_tri(N+1,rfq,sfq);
M = Vq'*diag(wq)*Vq;
Mf2 = Vfq'*diag(wfq)*Vfq2;

LIFT2 = M\Mf2(:,1:N+2);
LIFT2(abs(LIFT2)<1e-8) = 0;