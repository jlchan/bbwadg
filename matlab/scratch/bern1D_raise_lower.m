N = 5

% deriv ops

% degree raising
r = JacobiGL(0,0,N+1);
V1 = bern_basis_1D(N,r);
V2 = bern_basis_1D(N+1,r);
E1 = V2\V1; 
E1(abs(E1)<1e-8) = 0

% % undersample degree lowering
% r = JacobiGQ(0,0,N);
% V1 = bern_basis_1D(N,r);
% V2 = bern_basis_1D(N+1,r);
% E = V1\V2; 
% E(abs(E)<1e-8) = 0

% project down
[rq w] = JacobiGQ(0,0,N+1);
V1q = bern_basis_1D(N,rq);
V2q = bern_basis_1D(N+1,rq);
M2 = V2q'*diag(w)*V2q;
M1 = V1q'*diag(w)*V1q;
E2 = M1\(E1'*M2);
E2(abs(E2)<1e-8) = 0;

rp = linspace(-1,1,250);
V1p = bern_basis_1D(N,rp);
V2p = bern_basis_1D(N+1,rp);
c = randn(N+2,1);
plot(rp,V1p*E2*c,'b-');hold on;
plot(rp,V1p*E1'*c,'r-');
plot(rp,V2p*c,'k--');
rq = JacobiGQ(0,0,N);
hold on;plot(rq,rq*0,'o')

