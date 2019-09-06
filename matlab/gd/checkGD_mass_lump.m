% clear

N = 9;
K = 150;

% plotting points
Npts = 1000; a = -K/2; b = -a;
rp = linspace(a,b,Npts)';

[Vp VX] = GDVDM(N,K,rp);
VX = VX';
i = K/2+1; 

plot(VX,VX*0,'o')
hold on
plot(rp,sum(Vp(:,N+2:end-(N+2)),2),'linewidth',2);


re = linspace(-1,1,50)';
rK = [];
for e = 1:K
    h = (max(VX)-min(VX))/K;
    rK = [rK; h*(re+1)/2 + VX(e)];        
end

R = GDVDM(N,K,rK); phi = R(:,i);
p = 1+(rK).^2; %(rK-VX(i)).^(2);
moment = reshape(phi.*p,length(re),K);

mphi = sum(moment(:,end/2+1:end),2);
clf
% plot(re,mphi/max(abs(moment(:))),'o')
plot(re,moment,'o')
hold on
% plot(re,-moment(:,K/2+1),'o')

% 1st N+1 basis functions get modified - can orthogonalize w.r.t others?
% should just be able to modify first p+1 entries

%%

h = @(r) repmat(diff(VX),length(r),1);
map = @(r) reshape(h(r).*(repmat(r,1,K)+1)/2 + repmat(VX(1:end-1),length(r),1),length(r)*K,1);
L = (max(VX)-min(VX))/2;

% mass quadrature
rq = []; wq = [];
rqfull = []; wqfull = [];
for e = 1:K
    off = (N+1)/2-1; % this seems to work for N=1,3,5,7
    if e > off+1 && e < K - off % interior
%         [rqe wqe] = JacobiGL(0,0,1);
        [rqe wqe] = JacobiGQ(0,0,N);
    else
        [rqe wqe] = JacobiGQ(0,0,N);
    end
    h = (max(VX)-min(VX))/K;
    rq = [rq; h*(rqe+1)/2 + VX(e)];
    wq = [wq; h/2*wqe];
    
    % full quadrature
    [rqe wqe] = JacobiGQ(0,0,N);
    D1D = GradVandermonde1D(N,rqe)/Vandermonde1D(N,rqe);
    rqfull = [rqfull; h*(rqe+1)/2 + VX(e)];
    wqfull = [wqfull; h/2*wqe];    
end

Vqfull = GDVDM(N,K,rqfull);
Mfull = Vqfull'*diag(wqfull)*Vqfull;
VN = GDVDM(N,K,rqfull);
DN = kron(eye(K),D1D);
rx = 2;
Q = VN'*diag(wqfull)*(rx*DN*VN); % Galerkin first deriv mat
Q(abs(Q)<1e-8) = 0;

% possibly reduced quadrature versions
% [H, ~, X] = CSBPp4(K+1);
% wq = diag(H)*L; rq = X*L;
Vq = GDVDM(N,K,rq);
M = (Vq'*diag(wq)*Vq);
Pq = M\(Vq'*diag(wq));

[nm,~] = size(Mfull);
M2 = Mfull;
pskip = N+2;
inds = pskip:nm-(pskip)+1;
MDiag = diag(sum(M2,2));
M2(inds,:) = MDiag(inds,:);


f = @(x) exp(sin(pi*x/L));
df = @(x) (pi*cos((pi*x)/L).*exp(sin((pi*x)/L)))/L;

f = @(x) sin(pi*x/L);
df = @(x) cos(pi*x/L)*pi/L;

clf
% plot(rp,df(rp),'--')
hold on
plot(VX,(M\Q)*f(VX)-df(VX),'o','linewidth',2,'markersize',16)
% plot(VX,(diag(1./sum(M,2))*Q)*f(VX),'x--')
plot(VX,(M2\Q)*f(VX)-df(VX),'x','linewidth',2,'markersize',16)

iids = N+2:K+1-(N+1);
Msub = M(iids,iids);
[V D] = eig(M);
[lam p] = sort(diag(D),'descend');
V = V(:,p);
ii = N;

% bids = setdiff(1:size(M,2),iids);
% A = M(bids,bids);
% B = M(bids,iids);
% C = M(iids,iids);
% x = randn(size(M,2),1);
% u = x(bids);
% v = x(iids);
% % V(:,ii)'*M*V(:,ii) - V(:,ii)'*M2*V(:,ii)
% 
% return
% 
% [V D] = eig(abs(C));
% [~,id] = max(diag(D));
% v = VX*0;
% v(iids) = V(:,id)*sign(V(ceil(K/2),id));
% v = v/max(v);
% clf
% plot(VX,v,'o')
% hold on
% plot(VX(bids),0*VX(bids),'x')
% % plot(reshape(rq,N+1,K),reshape(DN*VN*v,N+1,K),'o-')
% plot(rp,Vp*v)
% LL = (K/2);
% plot(rp,(1+cos(pi*rp/LL))/2)
% 
% 
% return
% 
% uu = 1+0*exp(sin(pi*VX(:)/L));
% uI = uu; uI(bids) = 0;
% uB = uu; uB(iids) = 0;
% clf
% % plot(rp,Vp*uu)
% hold on
% plot(rp,(Vp*uI))
% plot(rp,(Vp*uB))
% plot(rp,(Vp*uB).*(Vp*uI)*4,'k.-')
% plot(VX,uB,'o')
% plot(rp,Vp(:,1:pskip-1))
% 
% BC = [B;C];
% 
% 
% return
% 
clf
Vf = GDVDM(N,K,[min(VX); max(VX)]);
B = diag([-1;1]);

Bperiodic = 0*Vf;
Bperiodic(1,1) = -1;
Bperiodic(1,end) = 1;
Bperiodic(2,1) = 1;
Bperiodic(2,end) = -1;
%.5*Vf'*B*
lam = eig(M2\(Q + .5*Vf'*B*Bperiodic));
plot(lam,'o')
return

clf
plot(rp,Vp*Pq*f(rq),'o')
hold on
plot(rp,f(rp))

