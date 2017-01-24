N = 10;
rp = linspace(-1,1,1000)';
Vp = Vandermonde1D(N,rp);
[r w] = JacobiGL(0,0,N);
V = Vandermonde1D(N,r);
Vr = GradVandermonde1D(N,r);
Dr = Vr/V;
M = inv(V*V');
% M = diag(w);
f = randn(N+1,1);
J = 1+randn(N+1,1);

B = diag([-1;zeros(N-1,1);1]);
rhs1 = diag(J)*M*(Dr*f);
rhs2 = diag(J)*B*f - (Dr*diag(J))'*M*f;

rhs1-rhs2


%% check 2D version
N = 4;
r1D = JacobiGL(0,0,N);
[r1Dq  w1Dq] = JacobiGQ(0,0,N);
[r s] = meshgrid(r1D);
r = r(:); s = s(:);

rfq = [r1Dq; ones(N+1,1); -r1Dq; -ones(N+1,1)];
sfq = [-ones(N+1,1); r1Dq; ones(N+1,1); -r1Dq];
wfq = repmat(w1Dq,4,1);

V = zeros((N+1)^2);
Vr = zeros((N+1)^2);
Vs = zeros((N+1)^2);
Vfq = zeros(length(wfq),(N+1)^2);
sk = 1;
for i = 0:N
    for j = 0:N
        V(:,sk) = JacobiP(r,0,0,i).*JacobiP(s,0,0,j);
        Vr(:,sk) = GradJacobiP(r,0,0,i).*JacobiP(s,0,0,j);
        Vs(:,sk) = JacobiP(r,0,0,i).*GradJacobiP(s,0,0,j);
        Vfq(:,sk) = JacobiP(rfq,0,0,i).*JacobiP(sfq,0,0,j);        
        sk = sk + 1;
    end
end
Dr = Vr/V;
Ds = Vs/V;
Vfq = Vfq/V;

M = inv(V*V'); 
% M = diag(sum(M,2)); % lump
Mf = Vfq'*diag(wfq)*Vfq;

f = randn(N+1); 
f = f(:);
f = (1+r).*(1+s).*(1-r).*(1-s);

rhs1 = M*Dr*f;
rhs2 = Mf*f - Dr'*M*f;

norm(rhs1-rhs2)