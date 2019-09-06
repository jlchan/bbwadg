clear
N = 15;
K = 64;

useSBP = 1;

alpha = 0; % 0 = upwinding

% plotting points
Npts = 100; a = -K/2; b = -a;
rp = linspace(a,b,Npts)';

[Vp VX] = GDVDM(N,K,rp);
VX = VX';

% mass quadrature
rq = []; wq = [];
h = 1; 
for e = 1:K    
    % full quadrature
    [rqe wqe] = JacobiGQ(0,0,N);
    D1D = GradVandermonde1D(N,rqe)/Vandermonde1D(N,rqe);
    rq = [rq; h*(rqe+1)/2 + VX(e)];
    wq = [wq; h/2*wqe];    
end

h = 2/K;

VN = GDVDM(N,K,rq);
DN = kron(eye(K),D1D);
rx = 2;
Q = VN'*diag(wq)*(rx*DN*VN); % Galerkin first deriv mat
Q(abs(Q)<1e-8) = 0;

% possibly reduced quadrature versions
Vq = GDVDM(N,K,rq);
M = (Vq'*diag(wq)*Vq);
M2 = M;
pskip = N+2;
inds = pskip:size(M,2)-(pskip)+1;
MDiag = diag(sum(M2,2));
M2(inds,:) = MDiag(inds,:);

b = ones(size(M,2),1);
A = M2\M;
b = M2\b;
x = b;
Lmax = max(eig(M2\M));
Lmin = min(eig(M2\M));
maxit = 50;
tol = 1e-8;

[x rvec rhist] = cheb_iter(A,b,x,Lmax,Lmin,tol,maxit);

x = b;
D = M2;
R = M-M2;    
for i = 1:maxit
    x = D\(b-R*x);
    rvecj(i)=norm(M*x-b);
end

semilogy(rvec,'o--')
hold on
semilogy(rvecj,'x--')

