% clear

N = 15;
K = 64;

% plotting points
Npts = 500;
a = -K/2;
b = -a;
rp = linspace(a,b,Npts)';

[Vp VX] = GDVDM(N,K,rp);
% VX = VX';
% plot(VX,VX*0,'o')
% hold on
% plot(rp,Vp(:,1:N+1),'linewidth',2);
% return

% 1st N+1 basis functions get modified - can orthogonalize w.r.t others?
% should just be able to modify first p+1 entries

h = @(r) repmat(diff(VX),length(r),1);
map = @(r) reshape(h(r).*(repmat(r,1,K)+1)/2 + repmat(VX(1:end-1),length(r),1),length(r)*K,1);
L = (max(VX)-min(VX))/2;

h = (max(VX)-min(VX))/K;

% mass quadrature
rq = []; wq = [];
for e = 1:K
    % full quadrature
    [rqe wqe] = JacobiGQ(0,0,N);
    D1D = GradVandermonde1D(N,rqe)/Vandermonde1D(N,rqe);
    rq = [rq; h*(rqe+1)/2 + VX(e)];
    wq = [wq; h/2*wqe];
end

Vq = GDVDM(N,K,rq);
M = (Vq'*diag(wq)*Vq);

VN = GDVDM(N,K,rq);
DN = kron(eye(K),D1D);
rx = 2;
Q = VN'*diag(wq)*(rx*DN*VN); % Galerkin first deriv mat
Q(abs(Q)<1e-8) = 0;

B = zeros(size(Q)); B(1,1) = -1; B(end,end) = 1;
Bout = zeros(K+1); Bout(1,K+1) = 1; Bout(K+1,1) = -1;
Bbc = Bout + B;
% Bbc = B; Bbc(1,end) = 1; Bbc(end,1) = -1;

Bm = zeros(K+1); Bm(1,1) = 1; Bm(K+1,K+1) = 1;
Bp = zeros(K+1); Bp(1,K+1) = 1; Bp(K+1,1) = 1;
a = 1;
Qbc = Q - .5*(Bbc - a*(Bp-Bm));

% CN1D = max(abs(eig(Qbc,M)));

[V D] = eig(M); [lam p] = sort(diag(D),'descend'); V = V(:,p);

P = eye(size(M))-V(:,K:K+1)*V(:,K:K+1)'; % project out large-norm modes

lam1D = eig((M\Qbc));
max(abs(lam1D))
% xlim([-.1 .1])

%%

% 2D experiments
Bjumpx = kron(M,Bp - Bm);
a = 1;
Q2D = kron(M,Q) - .5*kron(M,Bbc - a*(Bp-Bm));

M2D = kron(M,M);

[V2D D2D] = eig(Q2D,M2D);
[lam2d p]= sort(diag(D2D),'descend');
V2D = V2D(:,p);

[x y] = meshgrid(VX);

% mesh(x,y,reshape(real(V2D(:,1)),K+1,K+1))
CN2D = max(abs(lam2d))
% plot(lam2d,'o')

% acoustics
Z = zeros(size(M2D));
Qx = kron(M,Q) - .5*kron(M,Bbc);
Qy = kron(Q,M) - .5*kron(Bbc,M);

Bjumpy = kron(Bp - Bm,M);
Bjump = Bjumpx+Bjumpy;

a = 1/2;
Ah = [a*Bjump Qx Qy;
    Qx a*Bjumpx a*Bjumpy;
    Qy a*Bjumpx a*Bjumpy];

Mh = kron(eye(3),M2D);

lam2d_acous = eig(Ah,Mh);
CN2D_acous = max(abs(lam2d_acous))







