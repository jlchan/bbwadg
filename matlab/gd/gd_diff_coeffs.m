clear

p = 5;
K = 16;

% plotting points
rp = linspace(-1,1,50)';
[Vp VX] = GDVDM(p,K,rp);
VX = VX';
i = K/2+1;

[rq wq] = JacobiGQ(0,0,p);
rk = [];
wk = rk;
for e = 1:K
    h = (max(VX)-min(VX))/K;
    rk = [rk; h*(rq+1)/2 + VX(e)];
    wk = [wk; h/2*wq];
    D1D = GradVandermonde1D(p,rq)/Vandermonde1D(p,rq);
end
R = GDVDM(p,K,rk);

M = R'*diag(wk)*R;
DN = kron(eye(K),D1D);
Q = R'*diag(wk)*(2*DN*R); % Galerkin first deriv mat
Q(abs(Q)<1e-8) = 0;

S = (2*DN*R)'*diag(wk)*(2*DN*R); % Galerkin 2nd deriv

M2 = M;
[nm,~] = size(M);
pskip = p+2; inds = pskip:nm-(pskip)+1;
MDiag = diag(sum(M2,2));
M2(inds,:) = MDiag(inds,:);
% M = M2;

phi = R(:,i);
r = p;
f = rk.^r;
moment = reshape(phi.*f,length(rq),K);

norm(sum(moment,2))
% plot(rq,moment,'x--')
% hold on

re = linspace(-p,p,p+1)';
q = (p+1)/2;
cj = [];
sk = 1;
for j = -q:q-1
    cj(sk) = (-1)^j * nchoosek(p,q+j);
    sk = sk + 1;
end
cj = cj/factorial(p);

summ = 0;
prod = 1 + 0*rq;
for j = -q:q-1
    summ = summ + cj(j+q+1)*(rq+j+1).^(r-1);
    prod = prod.*(rq + j + 1);
end

[Hs, Qs, Xs] = CSBPp3(K+1);
Hs = Hs*K/2;
% x = VX;
% r = p+1;
% norm((Hs\Qs)*x.^r - r*x.^(r-1))

%% check dispersion error

mid = round(size(Q,2)/2);
aa = Q(mid,mid+(1:p)); %/ M2(mid,mid);

cj = Q(mid,mid+(-p:p));
dj = S(mid,mid+(-p:p));
cjs = Qs(mid,mid+(-p:p));
x = VX(mid+(-p:p));
x0 = VX(mid);

r = 6;
f = @(x) x.^(r);
df = @(x) r*x.^(r-1);
d2f = @(x) r*(r-1)*x.^(r-2);
% cj*f(x) - df(x0)
% cjs*f(x) - df(x0)
% dj*f(x)-d2f(x0)

aaSBP = Qs(mid,mid+(1:p));% / MSBP(mid,mid);

Npts = 100;
Lmin = 1/Npts;
Lmax = pi-1/Npts;
Lvec = linspace(Lmin,Lmax,Npts);
for i = 1:Npts
    ell = Lvec(i);
    ilamhist_interior(i) = 2*sum(aa.*sin(ell*(1:p)));
    ilamhist_interior_SBP(i) = 2*sum(aaSBP.*sin(ell*(1:p)));
end

figure(1)
loglog(Lvec,abs(ilamhist_interior-Lvec),'bo-','linewidth',2,'DisplayName','Lumped GD');
hold on
loglog(Lvec,abs(ilamhist_interior_SBP-Lvec),'rs-','linewidth',2,'DisplayName','SBP');
loglog(Lvec,Lvec.^(p+2),'k--','DisplayName','h^{p+2}')
loglog(Lvec,Lvec.^(2*p+1),'k--','DisplayName','h^{2p+1}')
title('Interior node dispersion')
legend show

