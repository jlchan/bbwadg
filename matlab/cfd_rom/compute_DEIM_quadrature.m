clear
load Usnap_smooth.mat

colids = 1:10:size(Usnap,2);
ids = 1:size(Usnap,1)/3;
offset = size(Usnap,1)/3;
Usnap = [Usnap(ids,colids) Usnap(ids + offset,colids) Usnap(ids + 2*offset,colids)];
Vsnap = [Vsnap(ids,colids) Vsnap(ids + offset,colids) Vsnap(ids + 2*offset,colids)];

VqK = kron(eye(K),Vq);
[U,S,~] = svd(VqK*[Usnap Vsnap],0);
% U = diag(1./sqrt(wJq(:)))*U;

Nmodes = 10;
[Uh,~,~] = svd([U(:,1:Nmodes) ones(size(Usnap,1),1)],0);
% plot(xq(:),Uh(:,1:Nmodes),'.-')
% return

%%

Nmodes = 5;
xq = linspace(-1,1,1000)';
wJq = [.5; ones(length(xq)-2,1); .5]; wJq = wJq/sum(wJq)*2;
Uh = Vandermonde1D(Nmodes-1,xq);

UDEIM(:,1) = Uh(:,1);
[~,id] = max(abs(Uh(:,1)));
p = id;
for j = 2:Nmodes
    r = Uh(:,j)-UDEIM*(UDEIM(p,:)\Uh(p,j));
    [~,id] = max(abs(r));
    p(j) = id;
    UDEIM = [UDEIM r];    
end
xDEIM = xq(p);
wDEIM = wJq(:)'*(UDEIM/UDEIM(p,:));
% wJq(:)'*(Uh/Uh(p,:))
return

plot(xDEIM,wDEIM,'o')
hold on
[xq wq] = JacobiGQ(0,0,Nmodes);
plot(xq,wq,'x')
% plot([-1,1],[0 0],'--')


% wDEIM = wJq(:)'*(UDEIM/UDEIM(p,:));

return
s = diag(S);
semilogy(s,'o')
hold on
semilogy(s(1:Nmodes),'rx')


% Ui = reshape(U(:,i),Nq*K,3);
% rho = reshape(Ui(:,1),N+1,K);
% m = reshape(Ui(:,2),N+1,K);
% E = reshape(Ui(:,3),N+1,K);
