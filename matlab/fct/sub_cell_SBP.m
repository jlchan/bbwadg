N = 3;
[r w] = JacobiGL(0,0,N);
rf = [-1;cumsum(w)-1];
V = Vandermonde1D(N,r);
Vr = GradVandermonde1D(N,r);
Dr = Vr/V;

Qr = diag(w)*Dr;
B = Qr+Qr';

e = ones(N+2,1);
DF = diag(e(2:end),1) - diag(e);
DF = DF(1:end-1,:);

IF = [DF;[-1 ones(1,N) 1]]\[Qr;sum(B,1)];
IF = IF + ones(N+2,1)*w'/sum(w);

u = exp(sin(pi*r));
% u = sin(pi*r);
% u = 1+r;

% plot(r,u,'o--')
% hold on
% plot(rf,IF*u,'x--')

% central averaging
IFV = .5*diag(ones(N+1,1))+.5*diag(ones(N,1),-1);
IFV(1,:) = 0; IFV(1,1) = 1; 
IFV(N+2,:) = 0; IFV(N+2,N+1) = 1; 


rp = linspace(-1,1,100)';
Vp = Vandermonde1D(N,rp)/V;
% plot(rp,Vp*Dr*u)
% hold on
% plot(rp,Vp*diag(1./w)*DF*IFV*u,'--')


QFV = DF*IFV;

e = null(Qr');
% e = 0*e;
e = ones(N+1,1);



% IN = pinv([QFV;ones(N+1,1)'])*[Qr;ones(N+1,1)'];
IN = pinv(QFV)*Qr + ones(N+2,1)*w'/sum(w);
plot(r,IN*u,'o--')
hold on
plot(r,u,'.--')

