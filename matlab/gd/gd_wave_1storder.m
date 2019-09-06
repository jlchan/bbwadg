clear

global D L tau nx 

N = 11;
K = 250;

% N = 11;
% K = 250;

CFL = .50;
FinalTime = 10.5;
massLump = 0;
tau = 1;

% plotting points
Npts = 100; a = -K/2; b = -a;
rp = linspace(a,b,Npts)';
[Vp VX] = GDVDM(N,K,rp);
VX = VX';

% mass quadrature
h = 1;
[rqe wqe] = JacobiGQ(0,0,N);
D1D = GradVandermonde1D(N,rqe)/Vandermonde1D(N,rqe);
rq = []; wq = [];
for e = 1:K
    rq = [rq; h*(rqe+1)/2 + VX(e)];
    wq = [wq; h/2*wqe];
end

h = 2/K;

VN = GDVDM(N,K,rq);
DN = kron(eye(K),D1D);
rx = 2;
Q = VN'*diag(wq)*(rx*DN*VN); % Galerkin first deriv mat
Q(abs(Q)<1e-8) = 0;

% possibly reduced quadrature versions
Vq = GDVDM(N,K,rq);
M = (Vq'*diag(wq)*Vq);

h = 1/K;
VX = (K/2+VX)*h;
rq = (K/2+rq)*h;
rp = (K/2+rp)*h;
M = h*M;

dt = CFL*h;
Nsteps = ceil(FinalTime/dt);
Nstages = N+2;
dt = FinalTime/Nsteps;
e = ones(K+1,1);

if massLump
    M2 = M;
    pskip = N+2;
    inds = pskip:size(M,2)-(pskip)+1;
    MDiag = diag(sum(M2,2));
    M2(inds,:) = MDiag(inds,:);
    M = M2; % lumping or not
    
end

if min(eig(M+M')/2)<1e-8 
    keyboard
end
% [M, Q, X] = CSBPp3(K+1);
% M = M/2;

D = M\Q;

B = zeros(K+1,2);
B(1,1) = 1;
B(K+1,2) = 1;
L = M\B;
nx = [-1;1];

D(abs(D)<1e-10) = 0;
D = sparse(D);

h*max(abs(eig(M\(Q))))

% %% advec eigs
% 
% A = zeros(K+1);
% U = zeros(K+1,1);
% for i = 1:(K+1)
%     U(i) = 1;    
%     rhsu = advecRHS(U);
%     U(i) = 0;
%     A(:,i) = rhsu(:);
% end
% lam = eig(A);
% plot(lam,'o')
% title(sprintf('max real part = %g\n',max(real(lam))))
% keyboard

%% wave eigs

A = zeros(2*(K+1));
U = zeros(K+1,2);
for i = 1:2*(K+1)
    U(i) = 1;
    p = U(:,1);
    u = U(:,2);
    U(i) = 0;
    [rhsp rhsu] = waveRHS(p,u);
    A(:,i) = [rhsp(:); rhsu(:)];
end
lam = eig(A);
plot(lam,'o')
title(sprintf('max real part = %g\n',max(real(lam))))
% xlim([-1,1]*1e-8)
keyboard

%% rk coeffs

m = 41;
u0 = @(x,t) sin(m/2*pi*x).*cos(m/2*pi*t);
% u0 = @(x) exp(-100*(x-.5).^2);

p = u0(VX,0);
u = zeros(K+1,1);

% plot(VX,(M\Q)*exp(VX)-exp(VX),'o')
% hold on
% plot(VX,(M2\Q)*exp(VX)-exp(VX),'x')
% plot(VX,exp(VX),'--')
% plot(VX,u0(VX))
% return

res = 0*u;
for i = 1:Nsteps   
    
    rhsp = p;
    rhsu = u;
    for j = 1:Nstages
        [rhsp rhsu] = waveRHS(rhsp,rhsu);
        rhsp = dt/j * rhsp;
        rhsu = dt/j * rhsu; 
        p = p + rhsp; 
        u = u + rhsu; 
    end
    
    if (mod(i,10)==0 || i==Nsteps)
        plot(VX,p-u0(VX,i*dt),'o--')
%         axis([0 1 -1e-4 1e-4])
        
%         plot(VX,p,'o--')
%         axis([0 1 -1.5 1.5])
        
        title(sprintf('time = %f, tstep %d/%d\n',dt*i,i,Nsteps))
        drawnow
    end
end

%norm(u-u0(mod(VX+FinalTime,1)),'inf')
% norm(p-u0(VX,FinalTime),'inf')
err2 = (Vq*p-u0(rq,FinalTime)).^2;
L2err = sqrt(sum(sum(err2)))

function rhsu = advecRHS(u)

global D L tau nx 

K = size(L,1)-1;

uM = u([1;K+1]);
uP = uM([2;1]);
ujump = (uP-uM);
f = .5*ujump.*nx - .5*tau*ujump;
rhsu = -((D*u) + L*f);
% rhsu = L*ujump;

end

function [rhsp rhsu] = waveRHS(p,u)

global D L tau nx 

K = size(L,1)-1;
pM = p([1;K+1]); 
uM = u([1;K+1]); 

% pP = pM([2;1]);
% uP = uM([2;1]);
pP = zeros(2,1);
uP = zeros(2,1);
pP(1) = -pM(1); % Dirichlet
pP(2) = pM(2);
uP(1) = uM(1);
uP(2) = -uM(2); % Neumann

pjump = (pP-pM);
ujump = (uP-uM);
fu = .5*pjump.*nx - .5*tau*(ujump.*nx).*nx;
fp = .5*ujump.*nx - .5*tau*pjump;
rhsp = -(D*u + L*fp);
rhsu = -(D*p + L*fu);

end


