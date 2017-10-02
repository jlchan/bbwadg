%% splines
clear

NB = 9;
Ksub = 1;
N = NB+Ksub-1;
smoothKnots = 0;

if Ksub==1
    [r1D w1D] = JacobiGL(0,0,N);
else
    VX = linspace(-1,1,Ksub+1);
    t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec
    for i = 1:N+1
        r1D(i,1) = mean(t((i+1):(i+NB))); % greville
    end
end
[BVDM M D1D R rq wq Vq] = bsplineVDM(NB,Ksub,r1D,smoothKnots); % VDM for interp, mass, M\S

wex = @(x) 1 + .5*sin(2*pi*x);
w = wex(rq);

M = Vq'*diag(wq)*Vq;
Mw = Vq'*diag(wq.*w)*Vq;
Mwadg = M*((Vq'*diag(wq./w)*Vq)\M);
[W D] = eig(Mw,Mwadg);
[lam,p] = sort(diag(D),'descend');
W = W(:,p);
max(lam)

% check cheb iter
b = Vq'*(wq.*sin(pi*rq));
A = Mw;
x = b;
P = eye(size(Mw));

P = inv(sqrtm(Mwadg));
A = P*A*P';
x = P*b;

% P = inv(Mwadg);
% A = P*A;
% x = P*b;

[x rvec] = cheb_iter(A,b,x,max(lam),min(lam),1e-11,25);
semilogy(rvec,'o--');
return
