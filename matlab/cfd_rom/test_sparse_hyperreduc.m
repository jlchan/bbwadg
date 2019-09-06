load Usnap_burgers_GLF.mat

K = size(Usnap,1);

xv = linspace(-1,1,K+1)';
x = .5*(xv(1:end-1)+xv(2:end));
dx = 2/K;

e = ones(K-1,1);
Q = diag(e,1)-diag(e,-1);
Q(1,:) = 0; Q(:,1) = 0;
Q(end,:) = 0; Q(:,end) = 0;
Q(1,end) = -1; Q(end,1) = 1;
Q(1,2) = 1; Q(2,1) = -1;
Q(K-1,K) = 1; Q(K,K-1) = -1;

[Vr Sr,~] = svd(Usnap,0);

N = 25;
sig = diag(Sr);
tol = sqrt(sum(sig(N+1:end).^2)/sum(sig.^2))
V = Vr(:,1:N);

[QV,SQV,~] = svd(Q*V,0);
QV = QV(:,diag(SQV)>1e-10);

Va = V;
Vb = V;
sk = 1;
Vmass = [];
for i = 1:size(Va,2)
    for j = 1:size(Vb,2)
        Vmass(:,sk) = Va(:,i).*Vb(:,j);
        sk = sk + 1;
    end
end
[Vmass,Smass,~] = svd(Vmass,0);
smass = diag(Smass);
smass_energy = sqrt(1 - (cumsum(smass.^2)./sum(smass.^2)));
Vmass = Vmass(:,smass_energy > tol);

w0 = dx*ones(size(x(:)));
[wr id] = get_empirical_cubature(Vmass,w0,tol,inf); % don't scale to reduce roundoff


plot(x(id),0*id,'o')
% [Vt,St,~] = svd([V ones(size(V,1),1) Q*V]);
% Vt = Vt(:,diag(St)>1e-12);
% Qt = Vt*Vt'*Q*Vt*Vt';
% u2 = Qt*g;
