clear
Globals1D;

N = 5;

K1D = 2;
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

StartUp1D;

R = zeros(N*K+1,(N+1)*K);
gids = 1:N+1;
for e = 1:K    
    for i = 1:N+1
        R(gids(i),i + (N+1)*(e-1)) = 1;
    end
    gids = gids + N;
end

re = linspace(-1,1,N+1)';
[rq wq] = JacobiGQ(0,0,N);
% Vq = bern_basis_1D(N,rq);
% Ve = Vandermonde1D(N,re); Vq = Vandermonde1D(N,rq)/Ve;
Vq = Vandermonde1D(N,rq)/V;
Mref = Vq'*diag(wq)*Vq;

MK = kron(spdiag(J(1,:)),Mref);
PK = kron(spdiag(J(1,:)),inv(Mref));
M = R*MK*R';
P = R*PK*R';

% cond(diag(1./sum(M,2))*M)

M1ref = [2 1; 1 2]/3;
Jx = diff(diag(1./sum(R,2))*R*x(:));

Kx = length(Jx);
ids = 1:2;
M1 = zeros(size(M));
for e = 1:Kx
    M1(ids,ids) = M1(ids,ids) + M1ref * Jx(e);
    ids = ids + 1;
end
clc
cond(M1\M)
cond(M)
cond(diag(1./sum(M,2))*M) % lumped
cond(diag(1./diag(M))*M)  % 

rp = linspace(-1,1,100);
Vp = Vandermonde1D(N,rp)/V;
xp = Vp*x;

invM = inv(M);

for i = 1:Np*K
    clf
    plot(xp,Vp*reshape(R'*invM(:,i),Np,K))
    pause
end

return
% M1ref

e = ones(N*K+1,1);
h = 1/(N*K);
M1 = 4*diag(e) + diag(e(2:end),1)+ diag(e(2:end),-1);
M1(1,1) = M1(1,1)/2; M1(end,end) = M1(end,end)/2;
M1 = h*M1/3;

cond(M1\M)


return

% cond(M1\M)
iM = R'*inv(M);

rp = linspace(-1,1,150)';
Vp = Vandermonde1D(N,rp)/V;
xp = Vp*x;
plot(xp,Vp*reshape(iM(:,5),N+1,K))