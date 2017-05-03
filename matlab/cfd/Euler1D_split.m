function Euler1D

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;

% Order of polymomials used for approximation
N = 7;
K1D = 16;
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

Nq = 2*N+2;
Nq = N;
[rq wq] = JacobiGQ(0,0,Nq);
Vq = Vandermonde1D(N,rq)/V;
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));
xq = Vq*x;

[rho m E] = vortexSolution(xq,0);
rho = Pq*rho;
m = Pq*m;
E = Pq*E;

time = 0;

% Runge-Kutta residual storage
resrho = zeros(Np,K);
resm = zeros(Np,K);
resE = zeros(Np,K);

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));

dt   = .1*xmin;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;

M = inv(V*V');
% LIFT = diag(1./sum(M,2))*M*LIFT; % SEM

wJq = diag(wq)*(Vq*J);

% outer time step loop

figure(1)
for tstep=1:Nsteps
    
    % low storage RK
    for INTRK = 1:5
        [rhsrho, rhsm, rhsE] = EulerRHS1D(rho, m ,E);
        
        resrho = rk4a(INTRK)*resrho + dt*rhsrho;
        resm   = rk4a(INTRK)*resm + dt*rhsm;
        resE   = rk4a(INTRK)*resE + dt*rhsE;
        
        rho = rho + rk4b(INTRK)*resrho;
        m   = m + rk4b(INTRK)*resm;
        E   = E + rk4b(INTRK)*resE;                
    end;
    
    % Increment time
    time = time+dt;
    
    if 1 && mod(tstep,5)==0
        
        clf
        hold on;
        plot(xp,Vp*rho,'-');      plot(x,rho,'o')
%         plot(xp,Vp*m  ,'-');      plot(x,m,'s')
%         plot(xp,Vp*E  ,'-');      plot(x,E,'^')
        hold off
        
%         axis([-.5 .5 -1 3])
        title(sprintf('Time = %f\n',time))
        drawnow
    end
    
    %enorm(tstep) = sum(sum(J.*(M*(rho.^2 + m.^2 + E.^2))));
    mq = Vq*m;
    rhoq = Vq*rho;
    Eq = Vq*E;
    uq = mq./rhoq;    
    
    u = mq./rhoq;
    p = (gamma-1.0)*(Eq - 0.5*rhoq.*u.^2);
    
    S(tstep) = sum(sum(wJq.*real(log(p.*rhoq.^(-gamma)))));
    
    
    if mod(tstep,1000)==0
        disp(sprintf('tstep = %d out of %d\n',tstep,Nsteps))
    end
end;

% [rhoq mq Eq] = vortexSolution(xq,FinalTime);
% err = wJq.*((Vq*rho-rhoq).^2 + (Vq*m-mq).^2 + (Vq*E-Eq).^2);
% err = sqrt(sum(err(:)))

figure(2);
plot(dt*(1:Nsteps),S)
hold on

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
    % % steady sine solution
    % rho = 1+.2*sin(pi*x); p = ones(size(rho));
    % u = 0*rho;
    % E = p/(gamma-1) + .5*rho.*u.^2;
    % m = rho.*u;
    
    % sine solution
    rho = 2 + sin(pi*(x - t));
    u = ones(size(x));
    p = ones(size(x));
    m = rho;
    E = p/(gamma-1) + .5*rho.*u.^2;

elseif opt==3    
    
    p = 1.*(x < 0) + .1.*(x > 0);
    rho = 1.*(x < 0) + .125.*(x > 0);
    u = zeros(size(x));
    m = rho.*u;
    E = p/(gamma-1) + .5*rho.*u.^2;
    
end



function [rhsrho, rhsm, rhsE] = EulerRHS1D(rho_in, m_in ,E_in)

% function [rhsrho, rhsm, rhsE] = EulerRHS1D(rho, m ,E)
% Purpose  : Evaluate RHS flux in 1D Euler

Globals1D;
global gamma Vq Pq

rho = Vq*rho_in;
m = Vq*m_in;
E = Vq*E_in;

% compute maximum velocity for LF flux
u = m./rho;
p = (gamma-1.0)*(E - 0.5*rho.*u.^2);
cvel = sqrt(gamma*p./rho); 
lm   = Pq*(abs(u) + cvel);

LFc   = reshape(max(lm(vmapP),lm(vmapM)),Nfp*Nfaces,K);

jump = @(u) reshape(u(vmapP)-u(vmapM),Nfp*Nfaces,K);
Dh = @(u) rx.*(Dr*u) + .5*LIFT*(Fscale.*jump(u).*nx);
    
opt = 1;
if opt==1
    
    rhof = Pq*m;
    mf   = Pq*(rho.*u.^2 + p);
    Ef   = Pq*((E+p).*u);

    rhsrho  = -Dh(rhof);
    rhsm    = -Dh(mf);
    rhsE    = -Dh(Ef);
    
    LFc = LFc*0;
    rhsrho = rhsrho + .5*LIFT*(Fscale.*LFc.*jump(rho_in));
    rhsm   = rhsm   + .5*LIFT*(Fscale.*LFc.*jump(m_in));
    rhsE   = rhsE   + .5*LIFT*(Fscale.*LFc.*jump(E_in));

elseif opt==2
    
    rhsrho  = -Dh(Pq*(rho.*u));
    rhsm    = -Dh(Pq*(p + .5*rho.*u.^2)) - .5*(Pq*(rho.*u.*(Vq*Dh(Pq*u))) + Pq*(u.*(Vq*Dh(Pq*(rho.*u)))));
    
    rhoue = E - .5*rho.*u.^2;
%     rhsE    = -Dh(Pq*((rhoue + .5*rho.*u.^2 + p).*u));
    rhsE  = Dh(Pq*((rhoue + p).*u)) + .5*Dh(Pq*(.5*rho.*u.^3)) ...
        + .5*Pq*(rho.*u.^2.*(Vq*Dh(Pq*(.5*u))) + u.*(Vq*Dh(Pq*(.5*rho.*u.^2)))); 
    rhsE = -rhsE;
    
%     rhsrho  = -Dh(rhof);
%     rhsm    = -Dh(mf);
%     rhsE    = -Dh(Ef);
    
    rhsrho = rhsrho + .5*LIFT*(Fscale.*LFc.*jump(rho_in));
    rhsm   = rhsm   + .5*LIFT*(Fscale.*LFc.*jump(m_in));
    rhsE   = rhsE   + .5*LIFT*(Fscale.*LFc.*jump(E_in));
    
elseif opt==3
    
end

return


