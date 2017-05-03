clear

N = 3;

r = JacobiGL(0,0,N);
rp = linspace(-1,1,50)';
[rq wq] = JacobiGQ(0,0,2*N+1);

wex = @(x) 1 + .5*sin(2*pi*x);
% wex = @(x) 1 + exp(3*sin(pi*x));
% wex = @(x) 1 + .5*x; 

w = wex(rq);

V = Vandermonde1D(N,r);
Vp = Vandermonde1D(N,rp);
Vq = Vandermonde1D(N,rq);
Pq = Vq'*diag(wq); 

Mw = Vq'*diag(wq.*w)*Vq;
invMwadg = Vq'*diag(wq./w)*Vq;
Mwadg = inv(invMwadg);
[W D] = eig(Mw,Mwadg);
[lam,p] = sort(diag(D),'descend');
W = W(:,p);
lam

min(Vp*Pq*w)
min(w)


vv = Vp*W;
vv = vv*diag(1./vv(1,:));
plot(rp,vv(:,:),'o')
hold on
W2 = Pq*(diag(1./w)*Vq);
vv = Vp*W2;
vv = vv*diag(1./vv(1,:));
plot(rp,vv(:,end))

%% splines
clear

NB = 3;
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
% wex = @(x) 1./log(2+x);
% wex = @(x) 1 + exp(3*sin(pi*x));
% wex = @(x) 1 + .5*x; 
w = wex(rq);

Mw = Vq'*diag(wq.*w)*Vq;
Mwadg = M*((Vq'*diag(wq./w)*Vq)\M);
[W D] = eig(Mw,Mwadg);
[lam,p] = sort(diag(D),'descend');
W = W(:,p);
max(lam)

%% check cheb iter

b = Vq'*(wq.*sin(pi*rq));
A = Mw;
x = b;
P = eye(size(Mw));
P = inv(Mw);
Lmax = max(lam);
Lmin = min(lam);

% cheb iter 
p = zeros(size(b));
d = (Lmax+Lmin)/2;
c = (Lmax-Lmin)/2;
beta = 0; alpha = 1.0; 
r = b-A*x;    
for i = 1:25
    z = P*r;
    p = r + beta*p;    
    
    alpha = 1/(d-beta/alpha);
    x = x + alpha*p;            
    beta = (c*alpha/2)^2;
    
    r = b-A*x;    

    rvec(i) = norm(r);      
end
semilogy(rvec,'o--')
%% compare GLL and exact mass

N = 4;
[r w] = JacobiGL(0,0,N);
V = Vandermonde1D(N,r);
M = inv(V*V');
MGLL = diag(w); % or diag(sum(M,2))
[W D]  = eig(M,MGLL)
lam = diag(D);


%% 2D
 
N = 3;
[r s] = Nodes2D(N); [r s] = xytors(r,s);
[rq sq wq] = Cubature2D(2*N+2);

wex = @(r,s) 1 + .75*sin(4*pi*r).*sin(4*pi*s);
w = wex(rq,sq);

V = Vandermonde2D(N,r,s);
Vq = Vandermonde2D(N,rq,sq)/V;
Pq = (V*V') * Vq'*diag(wq); 

Mw = Vq'*diag(wq.*w)*Vq;
invMwadg = (V*V')*Vq'*diag(wq./w)*Vq*(V*V');
[W D] = eig(invMwadg*Mw);
[lam,p] = sort(diag(D));
W = W(:,p);
[max(lam),min(lam)]

sort(eig(Mw,inv(invMwadg)))