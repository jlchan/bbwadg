clear
clc

N = 4;

r = JacobiGL(0,0,N);
rp = linspace(-1,1,75)';
[rq wq] = JacobiGQ(0,0,N+1);

wex = @(x) 1 + .75*sin(2*pi*x);
k = 1;
wex = @(x) 1 + exp(1+sin(k*pi*x));
wex = @(x) 1./(2 + x.^(N));

w = wex(rq);
f = exp(rq);

V = Vandermonde1D(N,r);
Vp = Vandermonde1D(N,rp)/V;
Vq = Vandermonde1D(N,rq)/V;
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq));

Mw = Vq'*diag(wq.*w)*Vq;
invMwadg = inv(M)*Vq'*diag(wq./w)*Vq*inv(M);
% Mwadg = inv(invMwadg);

uw = Mw\(Vq'*(wq.*f));
% u = (1:N+1)';
% uw = Mw\(M*u);


uwadg = invMwadg*(Vq'*(wq.*f));
uwadg0 = uwadg;
uwadg = uwadg + (wq'*(f - (w.*(Vq*uwadg)))/ (wq'*w));
% uwadg = invMwadg*(M*u);
% uwadg = uwadg - wq'*(w.*(Vq*uwadg))/(wq'*w) + ((wq'*(Vq*u))/ (wq'*w));

plot(rp,Vp*uw)
hold on
plot(rp,Vp*uwadg,'o--')


title(sprintf('diff in avg = %g, L2diff = %g\n',abs(wq'*(w.*(Vq*uw)) - wq'*(w.*(Vq*uwadg))),sqrt(wq'*((Vq*uwadg0-Vq*uwadg).^2))))



%%

clear
clc

N = 4;

r = JacobiGL(0,0,N);
rp = linspace(-1,1,75)';
[rq wq] = JacobiGQ(0,0,N+1);

k = 1;
w11 = @(x) 2 + exp(1+sin(k*pi*x));
w12 = @(x) 1 + x; 
w22 = @(x) 2 + exp(cos(k*pi*x));
% wex = @(x) 1 + 0*x; 

F = [exp(rq);exp(rq)];


V = Vandermonde1D(N,r);
Vp = Vandermonde1D(N,rp)/V;
Vq = Vandermonde1D(N,rq)/V;
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq));

W = @(x) [w11(x) w12(x);  w12(x) w22(x)];     
idetW = @(x) 1./(w11(x).*w22(x) -w12(x).^2);
invW = @(x) idetW(x) .* [w22(x) -w12(x);-w12(x) w11(x)];
iw11 = @(x) idetW(x) .* w22(x);
iw12 = @(x) idetW(x) .* -w12(x); 
iw22 = @(x) idetW(x) .* w11(x); 


Mw = [Vq'*diag(wq.*w11(rq))*Vq Vq'*diag(wq.*w12(rq))*Vq;
    Vq'*diag(wq.*w12(rq))*Vq Vq'*diag(wq.*w22(rq))*Vq];

Minvw = [Vq'*diag(wq.*iw11(rq))*Vq Vq'*diag(wq.*iw12(rq))*Vq;
    Vq'*diag(wq.*iw12(rq))*Vq Vq'*diag(wq.*iw22(rq))*Vq];
invMwadg = blkdiag(inv(M),inv(M))*Minvw*blkdiag(inv(M),inv(M));

Uw = Mw\(kron(eye(2),Vq'*diag(wq))*F);
Uwadg = invMwadg*(kron(eye(2),Vq'*diag(wq))*F);

% conservative correction
u1 = Vq*Uwadg(1:N+1);
u2 = Vq*Uwadg(N+2:end);

Wavg = [wq'*w11(rq) wq'*w12(rq);
    wq'*w12(rq) wq'*w22(rq)];

uosc = Wavg\(blkdiag(wq',wq')*(F - [w11(rq).*u1 + w12(rq).*u2;w12(rq).*u1 + w22(rq).*u2]));
Uwadg = Uwadg + [uosc(1)*ones(N+1,1);uosc(2)*ones(N+1,1)];

plot(rp,Vp*reshape(Uw,N+1,2))
hold on
plot(rp,Vp*reshape(Uwadg,N+1,2),'o--')

err = sqrt(sum(sum(diag(wq)*(Vq*(reshape(Uw-Uwadg,N+1,2)).^2))));
title(sprintf('L2 err = %g, conserv err = %g',err,sum(Mw*Uw)-sum(Mw*Uwadg)))
