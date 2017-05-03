function Euler1D

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;

% Order of polymomials used for approximation
N = 5;
K1D = 8;
FinalTime = 2;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

vmapP(1) = vmapM(end); % make periodic
vmapP(end) = vmapM(1);

rp = linspace(-1,1,25)';
Vp = Vandermonde1D(N,rp)/V;
xp=  Vp*x;

% Set initial conditions
global gamma Vq Pq
gamma = 1.4;

Nq = 2*N+1;
[rq wq] = JacobiGQ(0,0,Nq);
Vq = Vandermonde1D(N,rq)/V;
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));
xq = Vq*x;

%% splitting params

global dWdU11 dWdU12 dWdU13 dWdU22 dWdU23 dWdU33
global dfdW11 dfdW12 dfdW13 dfdW22 dfdW23 dfdW33

dWdU11 = @(rho,m,E)(1.0./(E.*rho.*2.0-m.^2).^2.*(gamma.*m.^4+m.^4+E.^2.*gamma.*rho.^2.*4.0-E.*gamma.*m.^2.*rho.*4.0))./rho;
dWdU12 = @(rho,m,E)m.^3.*1.0./(E.*rho.*2.0-m.^2).^2.*-2.0;
dWdU13 = @(rho,m,E)rho.*(E.*rho-m.^2).*1.0./(E.*rho.*2.0-m.^2).^2.*-4.0;
dWdU22 = @(rho,m,E)rho.*(E.*rho.*2.0+m.^2).*1.0./(E.*rho.*2.0-m.^2).^2.*2.0;
dWdU23 = @(rho,m,E)m.*rho.^2.*1.0./(E.*rho.*2.0-m.^2).^2.*-4.0;
dWdU33 = @(rho,m,E)rho.^3.*1.0./(E.*rho.*2.0-m.^2).^2.*4.0;

dfdW11 = @(W1,W2,W3)(2.0.^(1.0./(gamma+1.0)).*W2.*W3.*gamma.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*(-(W2.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0))-W1.*W3.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*2.0).*(gamma-1.0)).^(gamma./(gamma+1.0)))./((W2.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0))-W1.*W3.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*2.0).*(gamma-1.0));
dfdW12 = @(W1,W2,W3)-(((W2.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0))-W1.*W3.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*2.0).*(gamma-1.0).*(-1.0./2.0)).^(gamma./(gamma+1.0)).*(W2.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0))+W2.^2.*gamma.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0))-W1.*W3.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*2.0+W1.*W3.*gamma.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*2.0))./((W2.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0))-W1.*W3.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*2.0).*(gamma-1.0));
dfdW13 = @(W1,W2,W3)(2.0.^(1.0./(gamma+1.0)).*W1.*W2.*gamma.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*(-(W2.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0))-W1.*W3.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*2.0).*(gamma-1.0)).^(gamma./(gamma+1.0)))./((W2.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0))-W1.*W3.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*2.0).*(gamma-1.0));
dfdW22 = @(W1,W2,W3)(W2.*((W2.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0))-W1.*W3.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*2.0).*(gamma-1.0).*(-1.0./2.0)).^(gamma./(gamma+1.0)).*(W2.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*3.0-W2.^2.*gamma.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0))-W1.*W3.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*6.0+W1.*W3.*gamma.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*6.0))./(W3.*(W2.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0))-W1.*W3.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*2.0).*(gamma-1.0));
dfdW23 = @(W1,W2,W3)-(2.0.^(-(gamma.*2.0+1.0)./(gamma+1.0)).*1.0./W3.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*-2.0)./(gamma+1.0)).*(-(W2.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0))-W1.*W3.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*2.0).*(gamma-1.0)).^(gamma./(gamma+1.0)).*(W2.^4.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*4.0)./(gamma+1.0)).*3.0-W2.^4.*gamma.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*4.0)./(gamma+1.0)).*4.0+W2.^4.*gamma.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*4.0)./(gamma+1.0))+W1.^2.*W3.^2.*gamma.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*4.0)./(gamma+1.0)).*4.0-W1.*W2.^2.*W3.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*4.0)./(gamma+1.0)).*6.0-W1.^2.*W3.^2.*gamma.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*4.0)./(gamma+1.0)).*4.0-W1.*W2.^2.*W3.*gamma.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*4.0)./(gamma+1.0)).*4.0+W1.*W2.^2.*W3.*gamma.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*4.0)./(gamma+1.0)).*1.4e1))./((W2.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0))-W1.*W3.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*2.0).*(gamma-1.0));
dfdW33 = @(W1,W2,W3)(W2.*1.0./W3.^3.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*-2.0)./(gamma+1.0)).*((W2.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0))-W1.*W3.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*2.0).*(gamma-1.0).*(-1.0./2.0)).^(gamma./(gamma+1.0)).*(W2.^4.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*4.0)./(gamma+1.0))-W2.^4.*gamma.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*4.0)./(gamma+1.0)).*2.0+W2.^4.*gamma.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*4.0)./(gamma+1.0))+W1.^2.*W3.^2.*gamma.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*4.0)./(gamma+1.0)).*4.0-W1.*W2.^2.*W3.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*4.0)./(gamma+1.0)).*2.0-W1.^2.*W3.^2.*gamma.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*4.0)./(gamma+1.0)).*2.0-W1.*W2.^2.*W3.*gamma.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*4.0)./(gamma+1.0)).*4.0+W1.*W2.^2.*W3.*gamma.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*4.0)./(gamma+1.0)).*6.0))./((W2.^2.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0))-W1.*W3.*((W1.*W3-W2.^2.*(1.0./2.0)).*(gamma-1.0)).^(-(gamma+1.0)./(gamma-1.0)).^((gamma.*2.0)./(gamma+1.0)).*2.0).*(gamma-1.0));


%% 

[rhoq mq Eq] = vortexSolution(xq,0);
plot(xq,mq,'o')
[W1q W2q W3q] = UToW(rhoq,mq,Eq);
W1 = Pq*W1q;
W2 = Pq*W2q;
W3 = Pq*W3q;
% [rho m E] = WToU(W1q,W2q,W3q);
% hold on
% plot(xq,m,'x')
% return

time = 0;

% Runge-Kutta residual storage
res1 = zeros(Np,K);
res2 = zeros(Np,K);
res3 = zeros(Np,K);

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));

dt   = .1*xmin;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;

M = inv(V*V');
% LIFT = diag(1./sum(M,2))*M*LIFT; % SEM
global wJq
wJq = diag(wq)*(Vq*J);

% outer time step loop

figure(1)
for tstep=1:Nsteps
    
    ssprk = [1 .25 2/3];
    W1tmp = W1;
    W2tmp = W2;
    W3tmp = W3;
    for i = 1:3
        [rhs1, rhs2, rhs3] = EulerRHS1D(W1tmp,W2tmp,W3tmp);
        W1tmp = (1-ssprk(i))*W1 + ssprk(i)*(W1tmp + dt*rhs1);
        W2tmp = (1-ssprk(i))*W2 + ssprk(i)*(W2tmp + dt*rhs2);
        W3tmp = (1-ssprk(i))*W3 + ssprk(i)*(W3tmp + dt*rhs3);        
    end
    W1 = W1tmp;
    W2 = W2tmp;
    W3 = W3tmp;
    
%     % low storage RK
%     for INTRK = 1:5
%         [rhs1, rhs2, rhs3] = EulerRHS1D(W1,W2,W3);
%         
%         res1 = rk4a(INTRK)*res1 + dt*rhs1;
%         res2 = rk4a(INTRK)*res2 + dt*rhs2;
%         res3 = rk4a(INTRK)*res3 + dt*rhs3;
%         
% %         norm(res1)
% %         norm(res2)
% %         norm(res3)
%         W1 = W1 + rk4b(INTRK)*res1;
%         W2 = W2 + rk4b(INTRK)*res2;
%         W3 = W3 + rk4b(INTRK)*res3;                
%     end;
    
    % Increment time
    time = time+dt;
    
    if (1 && mod(tstep,10)==0) || (tstep==Nsteps)
        
        clf
        hold on;
        W1p = Vp*W1;
        W2p = Vp*W2;
        W3p = Vp*W3;
        [rhop mp Ep] = WToU(W1p,W2p,W3p);
        plot(xp,rhop,'-','linewidth',2);  
        hold off
        
        axis([-1 1 0 4])
%         axis([-.5 .5 -1 3])
        title(sprintf('Time = %f\n',time))
        drawnow
    end
    
    W1q = Vq*W1; W2q = Vq*W2; W3q = Vq*W3;
    [rho m E] = WToU(W1q,W2q,W3q);
    u = m./rho;
    p = (gamma-1.0)*(E - 0.5*rho.*u.^2);
    
    entropy = (p.*rho.^(-gamma)).^(1/(1+gamma));
    S(tstep) = sum(sum(wJq.*real(entropy)));
    
        
    if mod(tstep,1000)==0
        disp(sprintf('tstep = %d out of %d\n',tstep,Nsteps))
    end
end;

figure
semilogy(dt*(1:Nsteps),S)
keyboard

function [rho m E] = vortexSolution(x,t)

global gamma

opt = 1;
if opt==1
    rho = 2+exp(-(10*x).^2);
    % rho = 2+(x > 0);
    p = rho.^gamma;
    u = 0*rho;
    E = p/(gamma-1) + .5*rho.*u.^2;
    m = rho.*u;
    
elseif opt==2
    
    % steady sine solution
    rho = 1+.2*sin(pi*x); p = ones(size(rho));
    u = 0*rho;
    E = p/(gamma-1) + .5*rho.*u.^2;
    m = rho.*u;
    
%     % sine solution
%     rho = 2 + sin(pi*(x - t));
%     u = ones(size(x));
%     p = ones(size(x));
%     m = rho;
%     E = p/(gamma-1) + .5*rho.*u.^2;

elseif opt==3    
    
    p = 1.*(x < 0) + .1.*(x > 0);
    rho = 1.*(x < 0) + .125.*(x > 0);
    u = zeros(size(x));
    m = rho.*u;
    E = p/(gamma-1) + .5*rho.*u.^2;
    
end

function [W1 W2 W3] = UToW(rho,m,E)
global gamma 

u = m./rho;
p = (gamma-1.0)*(E - 0.5*rho.*u.^2);
pstar = -(rho.*p).^(-gamma/(gamma+1));

W1 = pstar.*E;
W2 = -pstar.*m;
W3 = pstar.*rho;


function [rho m E] = WToU(W1,W2,W3)
global gamma 

rhop = ((gamma-1)*(W1.*W3 - W2.^2/2)).^((1+gamma)/(1-gamma));
pstar = -(rhop).^(gamma/(gamma+1));
rho = pstar.*W3;
m = -pstar.*W2;
E = pstar.*W1;



function [rhs1, rhs2, rhs3] = EulerRHS1D(W1_in,W2_in,W3_in)

Globals1D;
global gamma Vq Pq

W1 = Vq*W1_in;
W2 = Vq*W2_in;
W3 = Vq*W3_in;

[rho m E] = WToU(W1,W2,W3);
u = m./rho;
p = (gamma-1.0)*(E - 0.5*rho.*u.^2);
cvel = sqrt(gamma*p./rho); 
lm   = Pq*(abs(u) + cvel);

LFc   = reshape(max(lm(vmapP),lm(vmapM)),Nfp*Nfaces,K);

jump = @(u) reshape(u(vmapP)-u(vmapM),Nfp*Nfaces,K);
Dh = @(u) rx.*(Dr*u) + .5*LIFT*(Fscale.*jump(u).*nx);

global dWdU11 dWdU12 dWdU13 dWdU22 dWdU23 dWdU33
global dfdW11 dfdW12 dfdW13 dfdW22 dfdW23 dfdW33

opt = 2;
if opt==1 % plain entropy variables
        
    dW1 = Vq*Dh(W1_in);
    dW2 = Vq*Dh(W2_in);
    dW3 = Vq*Dh(W3_in);
    
    % dF/dW * dW/dx
    %     r1 = dfdW11(rho,m,E).*dW1 + dfdW12(rho,m,E).*dW2 + dfdW13(rho,m,E).*dW3;
    %     r2 = dfdW12(rho,m,E).*dW1 + dfdW22(rho,m,E).*dW2 + dfdW23(rho,m,E).*dW3;
    %     r3 = dfdW13(rho,m,E).*dW1 + dfdW23(rho,m,E).*dW2 + dfdW33(rho,m,E).*dW3;
    r1 = dfdW11(W1,W2,W3).*dW1 + dfdW12(W1,W2,W3).*dW2 + dfdW13(W1,W2,W3).*dW3;
    r2 = dfdW12(W1,W2,W3).*dW1 + dfdW22(W1,W2,W3).*dW2 + dfdW23(W1,W2,W3).*dW3;
    r3 = dfdW13(W1,W2,W3).*dW1 + dfdW23(W1,W2,W3).*dW2 + dfdW33(W1,W2,W3).*dW3;

    % dW/dU * rhs
    rhs1 = dWdU11(rho,m,E).*r1 + dWdU12(rho,m,E).*r2 + dWdU13(rho,m,E).*r3;
    rhs2 = dWdU12(rho,m,E).*r1 + dWdU22(rho,m,E).*r2 + dWdU23(rho,m,E).*r3;
    rhs3 = dWdU13(rho,m,E).*r1 + dWdU23(rho,m,E).*r2 + dWdU33(rho,m,E).*r3;    
    
    rhs1 = -rhs1;
    rhs2 = -rhs2;
    rhs3 = -rhs3;
    
    LFc = LFc*0;
    rhs1 = Pq*rhs1 + .5*LIFT*(Fscale.*LFc.*jump(W1_in));
    rhs2 = Pq*rhs2 + .5*LIFT*(Fscale.*LFc.*jump(W2_in));
    rhs3 = Pq*rhs3 + .5*LIFT*(Fscale.*LFc.*jump(W3_in));

elseif opt==2
    
    F11 = dfdW11(W1,W2,W3);
    F12 = dfdW12(W1,W2,W3);
    F13 = dfdW13(W1,W2,W3);
    F22 = dfdW22(W1,W2,W3);
    F23 = dfdW23(W1,W2,W3);
    F33 = dfdW33(W1,W2,W3);
    
    A1 = Dh(Pq*(F11.*W1 + F12.*W2 + F13.*W3));
    A2 = Dh(Pq*(F12.*W1 + F22.*W2 + F23.*W3));
    A3 = Dh(Pq*(F13.*W1 + F23.*W2 + F33.*W3));
    
    dW1 = Vq*Dh(W1_in);
    dW2 = Vq*Dh(W2_in);
    dW3 = Vq*Dh(W3_in);
    B1 = Pq*(F11.*dW1 + F12.*dW2 + F13.*dW3);
    B2 = Pq*(F12.*dW1 + F22.*dW2 + F23.*dW3);
    B3 = Pq*(F13.*dW1 + F23.*dW2 + F33.*dW3);
    
    % split forms
    beta = (1+gamma)/(1-gamma);
    scale = 1 / (1+beta);    
    r1 = scale*Vq*(A1 + B1);
    r2 = scale*Vq*(A2 + B2);
    r3 = scale*Vq*(A3 + B3);
    
%     global wJq
%     sum(sum(wJq.*(W1.*r1 + W2.*r2 + W3.*r3)))

%     rhs1 = r1; rhs2 = r2; rhs3 = r3;

    % dW/dU * dF/dW * dW/dx
    rhs1 = dWdU11(rho,m,E).*r1 + dWdU12(rho,m,E).*r2 + dWdU13(rho,m,E).*r3;
    rhs2 = dWdU12(rho,m,E).*r1 + dWdU22(rho,m,E).*r2 + dWdU23(rho,m,E).*r3;
    rhs3 = dWdU13(rho,m,E).*r1 + dWdU23(rho,m,E).*r2 + dWdU33(rho,m,E).*r3;    
        
    LFc = LFc*0;
    rhs1 = Pq*rhs1 + .5*LIFT*(Fscale.*LFc.*jump(W1_in));
    rhs2 = Pq*rhs2 + .5*LIFT*(Fscale.*LFc.*jump(W2_in));
    rhs3 = Pq*rhs3 + .5*LIFT*(Fscale.*LFc.*jump(W3_in));
              
end

return




