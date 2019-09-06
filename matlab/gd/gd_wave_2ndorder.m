clear
% close all

N = 1;
K = 8;

massLump = 0;
computeEigs = 0;

CFL = .125;
FinalTime = 10.5;

% plotting points
Npts = 100; a = -K/2; b = -a;
[~, VX] = GDVDM(N,K,0);
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

M = h*M;
Ks = (1/h)*(rx*DN*VN)'*diag(wq)*(rx*DN*VN);
c2 = 1;
Ks = Ks*c2;

D = M\Q;
B = zeros(K+1,2);
B(1,1) = 1;
B(K+1,2) = 1;
L = M\B;
nx = [-1;1];

dt = CFL*h;
Nsteps = ceil(FinalTime/dt);
Nstages = N+2;
dt = FinalTime/Nsteps;
e = ones(K+1,1);

M2 = M;
pskip = N+2;
inds = pskip:size(M,2)-(pskip)+1;
MDiag = diag(sum(M2,2));
M2(inds,:) = MDiag(inds,:);
if massLump    
    M = M2; % lumping or not
end

if min(eig(M+M')/2)<1e-8
    keyboard
end

% lam = eig(M\(Ks+0*B*B'));
% [min(lam), max(lam)]
% return
lam = sort(eig(Ks,M),1,'ascend');
lamlump = sort(eig(Ks,M2),1,'ascend');
if ~isreal(lam) || ~isreal(lamlump)
    keyboard
end

fprintf('min full eig = %g, min lumped eig = %g\n',lam(1), lamlump(1))

%% initial cond

m = 41;
p0 = @(x,t) cos(m*pi*x).*cos(m*pi*t);
u0 = @(x,t) -m*pi*cos(m*pi*x).*sin(m*pi*t);

p0 = @(x,t) exp(-(x-.5).^2*400);

p = p0(VX,0);
u = u0(VX,0);
% plot(VX,p,'o--')
% hold on
% plot(VX,u,'x--')
% return

%% rk coeffs

rk4a = [            0.0 ...
    -567301805773.0/1357537059087.0 ...
    -2404267990393.0/2016746695238.0 ...
    -3550918686646.0/2091501179385.0  ...
    -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
    5161836677717.0/13612068292357.0 ...
    1720146321549.0/2090206949498.0  ...
    3134564353537.0/4481467310338.0  ...
    2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
    1432997174477.0/9575080441755.0 ...
    2526269341429.0/6820363962896.0 ...
    2006345519317.0/3224310063776.0 ...
    2802321613138.0/2924317926251.0];


resp = zeros(K+1,1);
resu = zeros(K+1,1);

a = 0;
b = a;

if computeEigs
    O = zeros(K+1);
    I = eye(K+1);
    II = eye(2*(K+1));
    MK = M\(Ks+B*diag([a;b])*B');
    A = [O MK;
        -I O];
    lam = eig(A);
    plot(lam,'o')
    
    mu = eig(MK);
    hold on
    plot([1i*sqrt(mu);-1i*sqrt(mu)],'x')
    
    title(sprintf('max real part = %5.5g\n',max(real(lam))))
    
    xlim([-1 1]*1e-4)
    return
end

D2 = M\Ks; D2(abs(D2)<1e-10) = 0; D2 = sparse(D2);

figure(1)
for i = 1:Nsteps
    
%     rhsp = p;
%     rhsu = u;
%     Nstages = N+1;
%     for j = 1:Nstages
% 
%         rhsu = -(D2*rhsp + L*([a;b].*rhsp([1;K+1])));
%         rhsp = rhsu;
%         
%         rhsu = dt/j * rhsu; 
%         rhsp = dt/j * rhsp;
%         
%         p = p + rhsp; 
%         u = u + rhsu; 
%     end
    
    for INTRK = 1:5
        rhsu = -(D2*p + L*([a;b].*p([1;K+1])));
        rhsp = u;
        
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resp = rk4a(INTRK)*resp + dt*rhsp;
        
        u = u + rk4b(INTRK)*resu;
        p = p + rk4b(INTRK)*resp;       
    end
    
    energy(i) = .5*(u'*M*u + p'*Ks*p);
    
    if (mod(i,50)==0 || i==Nsteps)
        %         plot(VX,p-u0(VX,i*dt),'o--')
        %         axis([0 1 -1e-4 1e-4])
        
        plot(VX,p,'o--')
        hold on
        plot(VX,p0(VX,i*dt),'x--')
        hold off
        axis([0 1 -1.5 1.5])
        
        title(sprintf('time = %f, tstep %d/%d\n',dt*i,i,Nsteps))
        drawnow
    end
end

figure(2)
semilogy(dt*(1:Nsteps),energy,'o--')
hold on

%norm(u-u0(mod(VX+FinalTime,1)),'inf')
% norm(p-u0(VX,FinalTime),'inf')
err2 = (Vq*p-p0(rq,FinalTime)).^2;
L2err = sqrt(sum(sum(err2)))

