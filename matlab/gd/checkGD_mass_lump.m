clear

N = 3;
K = 8;

% plotting points
Npts = 1000;
a = -K/2;
b = -a;
rp = linspace(a,b,Npts)';

[Vp VX] = GDVDM(N,K,rp);
VX = VX';
plot(VX,VX*0,'o')
hold on
plot(rp,Vp);return

%%

h = @(r) repmat(diff(VX),length(r),1);
map = @(r) reshape(h(r).*(repmat(r,1,K)+1)/2 + repmat(VX(1:end-1),length(r),1),length(r)*K,1);
L = (max(VX)-min(VX))/2;

if 1
    % quadrature
    [rqe wqe] = JacobiGQ(0,0,N);
    rq = [];% zeros(length(rqe),K);
    wq = rq;
    for e = 1:K
        if e==1 || e==K
            [rqe wqe] = JacobiGL(0,0,N);
        else
            [rqe wqe] = JacobiGL(0,0,N);
        end
        h = (max(VX)-min(VX))/K;
        rq = [rq; h*(rqe+1)/2 + VX(e)];
        wq = [wq; h/2*wqe];
        D1D{e} = GradVandermonde1D(N,rqe)/Vandermonde1D(N,rqe);
    end
    % rq = rq(:); wq = wq(:);
    % rq = map(rqe); rq = rq(:);
    % wq = repmat(wqe,1,K).*h(rqe)/2; wq = wq(:);
end

[H2 Q2 X2] = CSBPp2(K+1);
[H3 Q3 X3] = CSBPp3(K+1);
[H4 Q4 X4] = CSBPp4(K+1);
% wq = diag(H);
% rq = X*L;
% plot(rq,rq*0,'o');return

Vq = GDVDM(N,K,rq);
M = (Vq'*diag(wq)*Vq);
Pq = M\(Vq'*diag(wq));

VN = GDVDM(N,K,rq);
DN = blkdiag(D1D{:});%kron(eye(K),GradVandermonde1D(N,rqe)/Vandermonde1D(N,rqe)); % block diff
rx = 2;
Q = VN'*diag(wq)*(rx*DN*VN); % Galerkin first deriv mat
Q(abs(Q)<1e-8) = 0;

f = @(x) exp(sin(pi*x/L));
df = @(x) (pi*cos((pi*x)/L).*exp(sin((pi*x)/L)))/L;

f = @(x) x;
df = @(x) 1+0*x;

% plot(rq,(rx*DN*VN)*Pq*f(rq),'.')
%plot(rq,Vq*(M\Q)*Pq*f(rq),'.')
plot(VX,(M\Q)*Pq*f(rq),'.')
hold on;plot(VX,(diag(1./sum(M,2))*Q)*f(X2),'o')
% plot(rp,Vp*((1./diag(M)).*(Q*f(VX(:)))),'.')
hold on
plot(rp,df(rp),'-')
return

% check mass lumping
if 0
    [nm,~] = size(M);
    MDiag  = diag(sum(M,2));
    pskip = N+2;
    inds = pskip:nm-(pskip)+1;
    M(inds,:) = MDiag(inds,:);
end

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