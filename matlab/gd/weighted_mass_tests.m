N = 11;
K = 150;

Npts = 100;
rp = linspace(-K/2,K/2,Npts)';
[Vp, VX] = GDVDM(N,K,rp);  VX = VX(:);
L = max(VX)-min(VX);
h = VX(2)-VX(1);

% mass quadrature
[rqe wqe] = JacobiGQ(0,0,N);
rq = []; wq = [];
for e = 1:K
    rq = [rq; h*(rqe+1)/2 + VX(e)];
    wq = [wq; h/2*wqe];
end
Vq = GDVDM(N,K,rq);
M = (Vq'*diag(wq)*Vq);

M2w = M;
pskip = N+2;
inds = pskip:size(M,2)-(pskip)+1;
MDiag = diag(sum(M,2));
M2w(inds,:) = MDiag(inds,:);
% M = M2;

iids = N+2:K+1-(N+1);
Msub = M(iids,iids);
[V D] = eig(M);
[lam p] = sort(diag(D),'descend');
V = V(:,p);
ii = N;

bids = setdiff(1:size(M,2),iids);
A = M(bids,bids);
B = M(bids,iids);
C = M(iids,iids);
x = randn(size(M,2),1);
u = x(bids);
v = x(iids);
Mp = [A B;
    B' C];

%% check weighted version

w = @(x) 1+0*.25*sin(2*pi*x/L); % 1e-12 + rand(size(x));
wwqVX = w([rq;VX]);
wwq = wwqVX(1:length(rq));
wVX = wwqVX(length(rq)+1:end);

M = Vq'*diag(wq)*Vq;
Mw = Vq'*diag(wq.*wwq)*Vq;
Minvw = Vq'*diag(wq./wwq)*Vq;

f = @(x) exp(sin(pi*x/L));
bw = Vq'*diag(wq.*wwq)*f(rq);

% exact weighted mass inverse
u = Mw\bw;

% wadg
uwadg = M\(Minvw*(M\bw));

% lumped mass
M2w = Mw;
pskip = N+2;
inds = pskip:size(M,2)-(pskip)+1;
% MDiag = diag(sum(Mw,2));
MDiag = diag(wVX);
M2w(inds,:) = MDiag(inds,:);
ulump = M2w\bw;

% wadg for boundary block
A = Mw(bids,bids);
B = Mw(bids,iids);
C = Mw(iids,iids);
MA = Vq(:,bids)'*diag(wq)*Vq(:,bids);
MinvwA = Vq(:,bids)'*diag(wq./wwq)*Vq(:,bids);
MwadgA = MA*(MinvwA\MA);
A = MwadgA;
D = diag(w(VX(iids)));
% D = diag(sum(C,2));
invMw = [inv(A) -(A\B*inv(D)); 0*B' inv(D)];
invMw([bids,iids],[bids,iids]) = invMw;
ulumpwadg = invMw*bw;

% f = @(x) 0*x;
% hold on
semilogy(rp,abs(f(rp)-Vp*u),'bo--')
hold on
plot(rp,abs(f(rp)-Vp*uwadg),'rx-')
plot(rp,abs(f(rp)-Vp*ulump),'m^--')
plot(rp,abs(f(rp)-Vp*ulumpwadg),'k.-')
legend('Weighted','WADG','Lumped','Lumped WADG')
return

%% check for random weights

clear

while 1
    N = 3;
    K = 16;
    
    Npts = 500;
    rp = linspace(-K/2,K/2,Npts)';
    [Vp, VX] = GDVDM(N,K,rp);  VX = VX(:);
    L = max(VX)-min(VX);
    h = VX(2)-VX(1);
    
    % mass quadrature
    [rqe wqe] = JacobiGQ(0,0,N);
    rq = []; wq = [];
    for e = 1:K
        rq = [rq; h*(rqe+1)/2 + VX(e)];
        wq = [wq; h/2*wqe];
    end
    Vq = GDVDM(N,K,rq);
    M = (Vq'*diag(wq)*Vq);
    
    M2w = M;
    pskip = N+2;
    inds = pskip:size(M,2)-(pskip)+1;
    MDiag = diag(sum(M,2));
    M2w(inds,:) = MDiag(inds,:);
    
    
    w = @(x) 1e-12 + rand(size(x));
    wwqVX = w([rq;VX]);
    wwq = wwqVX(1:length(rq));
    wVX = wwqVX(length(rq)+1:end);
    
    Mw = Vq'*diag(wq.*wwq)*Vq;
    Minvw = Vq'*diag(wq./wwq)*Vq;
    
    % W = diag(wVX); invW = diag(1./wVX);
    % invW = diag(1./sum(Mw,2));
    %     invW = diag(sum(Minvw,2));
    %     W = diag(sum(Mw,2));
    
    %     M2w = M2*inv(Minvw)*M2;
    %     M2w = M*inv(Minvw)*M;
    M2w = Mw;
    MDiag = diag(sum(M2w,2));
    M2w(inds,:) = MDiag(inds,:);
    
    %     M2w(inds,:) = W(inds,:);
    
    lam = (eig(.5*(M2w + M2w')));
    fprintf('min eig = %f, max eig = %f\n',min(lam), max(lam))
    
    if min(real(lam))<0
        [V D] = eig(.5*(M2w + M2w'));
        x = V(:,find(diag(D)<0,1,'first'));
        Mw = M*inv(Minvw)*M;
        p = [setdiff(1:size(M2w,2),inds) inds];
        Mw = Mw(p,p);
        M2w = M2w(p,p);
        x = x(p);
        u = x(1:end-1); v = x(end); % assumes inds = size 1
        keyboard
    end
    
    %     u = exp(sin(2*pi*VX/L));
    %     plot(rp,Vp*u)%,'linewidth',2)
    %     hold on
    %     plot(rp,Vp*(M2w\(Mw*u)),'--')%,'linewidth',2)
end

