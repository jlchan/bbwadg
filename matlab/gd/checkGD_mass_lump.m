% clear

N = 3;
K = 64;

% plotting points
Npts = 1000;
a = -K/2;
b = -a;
rp = linspace(a,b,Npts)';

[Vp VX] = GDVDM(N,K,rp);
VX = VX';
plot(VX,VX*0,'o')
hold on
plot(rp,Vp(:,N+2),'linewidth',2);
% return

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
    if e==off+1 || e==K-off
        [rqe wqe] = JacobiGQ(0,0,N);
    elseif e > off+1 && e < K - off
        [rqe wqe] = JacobiGL(0,0,1);
%         [rqe wqe] = JacobiGQ(0,0,N);
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

% [H D tL tR X] = HGTpEQ3(K);
% [H Q X] = CSBPp2(K+1);
% [H3 Q3 X3] = CSBPp3(K+1);
% [H Q X] = CSBPp4(K+1);
% wq = diag(H); rq = X*L;
% plot(rq,rq*0,'o');return

Vq = GDVDM(N,K,rq);
M = (Vq'*diag(wq)*Vq);
Pq = M\(Vq'*diag(wq));

Vqfull = GDVDM(N,K,rqfull);
Mfull = Vqfull'*diag(wqfull)*Vqfull;

VN = GDVDM(N,K,rqfull);
DN = kron(eye(K),D1D);
rx = 2;
Q = VN'*diag(wqfull)*(rx*DN*VN); % Galerkin first deriv mat
Q(abs(Q)<1e-8) = 0;

[nm,~] = size(M);
M2 = Mfull;
pskip = N+2;
inds = pskip:nm-(pskip)+1;
MDiag = diag(sum(M,2));
M2(inds,:) = MDiag(inds,:);
% M2 = M2';

% M = M-diag(sum(M,2))+diag(sum(Mfull,2));

% offd = sum(M2(1:N+1,pskip:end),2);
% M2(1:N+1,pskip:end) = 0;
% M2(1:N+1,1:N+1) = M2(1:N+1,1:N+1) + diag(offd);
% offd = sum(M2(end-N:end,1:end-N-1),2);
% M2(end-N:end,1:end-N-1) = 0;
% M2(end-N:end,end-N:end) = M2(end-N:end,end-N:end) + diag(offd);
% clf;imagesc(M2);return

f = @(x) exp(sin(pi*x/L));
df = @(x) (pi*cos((pi*x)/L).*exp(sin((pi*x)/L)))/L;

% f = @(x) 1 + x;
% df = @(x) 1 + 0*x;

clf
% plot(rp,df(rp),'--')
hold on
plot(VX,(M\Q)*f(VX)-df(VX),'o','linewidth',2,'markersize',16)
% plot(VX,(diag(1./sum(M,2))*Q)*f(VX),'x--')
plot(VX,(M2\Q)*f(VX)-df(VX),'x','linewidth',2,'markersize',16)
return



clf
Vf = GDVDM(N,K,[min(VX); max(VX)]);
B = diag([-1;1]);

Bperiodic = 0*Vf;
Bperiodic(1,1) = -1;
Bperiodic(1,end) = 1;
Bperiodic(2,1) = 1;
Bperiodic(2,end) = -1;
%.5*Vf'*B*
lam = eig(M\(Q + .5*Vf'*B*Bperiodic));
plot(lam,'o')
return

clf
plot(rp,Vp*Pq*f(rq),'o')
hold on
plot(rp,f(rp))

function [V VX] = GDVDM(N,K,r)

nGhost = (N-1)/2;
Ka = (K/2 + nGhost); % half span + ghost point
VX = -Ka:Ka;

Npts = length(r);

% build ghost VDM
Np = length(-K/2:K/2);
V = zeros(Npts,Np);
for j = VX
    jid = Ka+1+j;
    for i = 1:Npts
        V(i,jid) = phi(r(i)-VX(jid),N);
    end
end

% modify extrapolation
ec = getExtrapCoeffs( N+1 );
for ig = nGhost:-1:1
    ileft  = 1+nGhost;
    iGhost = ileft-ig;
    for k = 2:length(ec)
        V(:,iGhost+k-1) = V(:,iGhost+k-1) + ec(k)*V(:,iGhost);
    end
end
for ig = nGhost:-1:1
    iright = K + nGhost + 1;
    iGhost = iright+ig;
    for k = 2:length(ec)
        V(:,iGhost-k+1) = V(:,iGhost-k+1)+ec(k)*V(:,iGhost);
    end
end
V = V(:,nGhost + (1:Np)); % extract non-ghost columns
VX = VX(nGhost+(1:Np));

end


function z = phi(xi,p)

z = 0;
a = -(p+1)/2;
b = -a;
if( xi >= a & xi < b )
    ia = p-ceil(xi-a)+1;
    ib = -(p-ceil(b-xi)+1);
    lp = [ia:-1:1,-1:-1:ib];
    
    den = 1;
    num = 1;
    for k = 1:length(lp)
        num = num*(xi+lp(k));
        den = den*lp(k);
    end
    z = num/den;
end

end