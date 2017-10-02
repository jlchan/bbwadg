function Euler1D

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;

% Order of polymomials used for approximation
N = 7;
K1D = 4;
FinalTime = 50;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(0,10,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

vmapP(1) = vmapM(end); % make periodic
vmapP(end) = vmapM(1);

rp = linspace(-1,1,25)';
Vp = Vandermonde1D(N,rp)/V;
xp=  Vp*x;

global VB DB VBe W
re = linspace(-1,1,N+1)';
[VB VBr] = bern_basis_1D(N,r);
DB = VB\VBr;
DB(abs(DB)<1e-8) = 0;
VBe = bern_basis_1D(N,re);
xe = Vandermonde1D(N,re)/V * x;

global Vq Pq
[rq wq] = JacobiGQ(0,0,N+1);
Vq = Vandermonde1D(N,rq)/V;
Pq = V*V' * Vq'*diag(wq);
xq = Vq*x;

% Set initial conditions
global gamma
gamma = 1.4;

[rhoq mq Eq uq pq] = vortexSolution(xq,0);
rho = Pq*rhoq;
m = Pq*mq;
E = Pq*Eq;
u = Pq*uq;
p = Pq*pq;

% for time = 0:.1:1;
%     [rho m E] = vortexSolution(xp,time);
%     plot(xp,rho)
%     hold on
% %     plot(xp,m)
% %     plot(xp,E)
%     hold off
%     drawnow
% end
% return

time = 0;

% Runge-Kutta residual storage
resrho = zeros(Np,K);
resm = zeros(Np,K);
resE = zeros(Np,K);

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));

dt   = .25*xmin;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;

M = inv(V*V');
wJq = diag(wq)*(Vq*J);

% outer time step loop
for tstep=1:Nsteps
    
    % low storage RK
    for INTRK = 1:5
        
%         [rhsrho, rhsm, rhsE] = EulerRHS1D_prim(rho, u ,p);                
        [rhsrho, rhsm, rhsE] = EulerRHS1D(rho, m ,E);        
        resrho = rk4a(INTRK)*resrho + dt*rhsrho;
        resm   = rk4a(INTRK)*resm + dt*rhsm;
        resE   = rk4a(INTRK)*resE + dt*rhsE;        
        rho = rho + rk4b(INTRK)*resrho;
        m   = m + rk4b(INTRK)*resm;
        E   = E + rk4b(INTRK)*resE;                        

% rho = rho + rk4b(INTRK)*resrho;
% u   = u + rk4b(INTRK)*resm;
% p   = p + rk4b(INTRK)*resE;
    end;
    
    % Increment time
    time = time+dt;
    
    if 1 && mod(tstep,5)==0
        
        clf
        hold on;
        plot(xp,Vp*rho,'-')
        plot(x,rho,'o')
        plot(xp,Vp*m,'-')
        plot(x,m,'o')
        plot(xp,Vp*E,'-')
        plot(x,E,'o')
        hold off
        
        ylim([-1 7])
        title(sprintf('Time = %f\n',time))
        drawnow                
    end
        
    rhoq = Vq*rho;
    if min(rhoq(:))<1e-8
        keyboard
    end
    e2 = (Vq*m).^2./rhoq;    
    enorm(tstep) = sqrt(wJq(:)'*e2(:)); % kinetic energy
    
    if mod(tstep,1000)==0
        disp(sprintf('tstep = %d out of %d\n',tstep,Nsteps))
    end
end;

figure;
plot(enorm)
keyboard


function [rho m E u p] = vortexSolution(x,t)

global gamma

x0 = 5; 
rho = 5 + exp(-(x-x0).^2); 
% rho = 1 + (x > x0);
p = rho.^gamma;
u = 0*rho;
E = p/(gamma-1) + .5*rho.*u.^2;
m = rho.*u;

% % steady sine solution
% rho = 1+.2*sin(pi*x); p = ones(size(rho)); 
% u = 0*rho;
% E = p/(gamma-1) + .5*rho.*u.^2;
% m = rho.*u; 

% sine solution
rho = 2 + sin(pi*x);
u = ones(size(x));
p = ones(size(x));
m = rho;
E = p/(gamma-1) + .5*rho.*u.^2;



function [rhsrho, rhsm, rhsE] = EulerRHS1D(rho, m ,E)

% function [rhsrho, rhsm, rhsE] = EulerRHS1D(rho, m ,E)
% Purpose  : Evaluate RHS flux in 1D Euler

Globals1D;
global gamma
global Vq Pq

% compute maximum velocity for LF flux
u = Pq*((Vq*m)./(Vq*rho));
rhou2 = Pq*((Vq*rho).*(Vq*u).^2);
pres = (gamma-1.0)*(E - 0.5*rhou2);
cvel = sqrt(gamma*Pq*((Vq*pres)./(Vq*rho))); 
lm   = abs(u) + cvel;
lm = zeros(size(u)); 

% Compute fluxes
rhof = m; 
% mf   = rho.*u.^2 + pres; 
mf = rhou2 + pres; 
%Ef   = (E+pres).*u;
Ef   = Pq*(Vq*(E+pres).*(Vq*u));

% Compute jumps at internal faces
drho = zeros(Nfp*Nfaces,K);  drhof = zeros(Nfp*Nfaces,K); 
dm   = zeros(Nfp*Nfaces,K);  dmf   = zeros(Nfp*Nfaces,K); 
dEf  = zeros(Nfp*Nfaces,K);  dE    = zeros(Nfp*Nfaces,K); 

drho(:)  = rho(vmapM)-  rho(vmapP);
dm(:)    = m(vmapM)- m(vmapP);
dE(:)    = E(vmapM)- E(vmapP);
drhof(:) = rhof(vmapM)- rhof(vmapP);
dmf(:)   = mf(vmapM)-mf(vmapP);
dEf(:)   = Ef(vmapM)-Ef(vmapP);

LFc    = zeros(Nfp*Nfaces,K); 
LFc(:) = max(lm(vmapP),lm(vmapM));

% Compute fluxes at interfaces
drhof(:) = nx(:).*drhof(:) - LFc(:).*drho(:);
dmf(:)   = nx(:).*dmf(:)   - LFc(:).*dm(:);
dEf(:)   = nx(:).*dEf(:)   - LFc(:).*dE(:);

% compute right hand sides of the PDE?s
rhsrho  = -rx.*(Dr*rhof) + LIFT*(Fscale.*drhof/2);
rhsm    = -rx.*(Dr*mf)   + LIFT*(Fscale.*dmf/2);
rhsE    = -rx.*(Dr*Ef)   + LIFT*(Fscale.*dEf/2);

return


% function [rhsrho, rhsm, rhsE] = EulerRHS1D_prim(rho, u ,p)
% 
% % function [rhsrho, rhsm, rhsE] = EulerRHS1D(rho, m ,E)
% % Purpose  : Evaluate RHS flux in 1D Euler
% 
% Globals1D;
% global gamma
% global Vq Pq
% 
% % compute maximum velocity for LF flux
% u = Pq*((Vq*m)./(Vq*rho));
% rhou2 = Pq*((Vq*rho).*(Vq*u).^2);
% E = p/(gamma-1) + .5*rhou2;
% cvel = sqrt(gamma*Pq*((Vq*p)./(Vq*rho))); 
% lm   = abs(u) + cvel;
% lm = zeros(size(u)); 
% 
% % Compute fluxes
% rhof = Pq*((Vq*rho).*(Vq*u));
% mf = rhou2 + p; 
% %Ef   = Pq*(Vq*(E+pres).*(Vq*u));
% Ef   = Pq*(Vq*(E+p).*(Vq*u));
% 
% % Compute jumps at internal faces
% drho = zeros(Nfp*Nfaces,K);  drhof = zeros(Nfp*Nfaces,K); 
% dm   = zeros(Nfp*Nfaces,K);  dmf   = zeros(Nfp*Nfaces,K); 
% dEf  = zeros(Nfp*Nfaces,K);  dE    = zeros(Nfp*Nfaces,K); 
% 
% drho(:)  = rho(vmapM)-  rho(vmapP);
% dm(:)    = m(vmapM)- m(vmapP);
% dE(:)    = E(vmapM)- E(vmapP);
% drhof(:) = rhof(vmapM)- rhof(vmapP);
% dmf(:)   = mf(vmapM)-mf(vmapP);
% dEf(:)   = Ef(vmapM)-Ef(vmapP);
% 
% LFc    = zeros(Nfp*Nfaces,K); 
% LFc(:) = max(lm(vmapP),lm(vmapM));
% 
% % Compute fluxes at interfaces
% drhof(:) = nx(:).*drhof(:) - LFc(:).*drho(:);
% dmf(:)   = nx(:).*dmf(:)   - LFc(:).*dm(:);
% dEf(:)   = nx(:).*dEf(:)   - LFc(:).*dE(:);
% 
% % compute right hand sides of the PDE?s
% rhsrho  = -rx.*(Dr*rhof) + LIFT*(Fscale.*drhof/2); % d(rho)/dt
% rhsm    = -rx.*(Dr*mf)   + LIFT*(Fscale.*dmf/2); % d(rho*u)/dt = rho * d(u)/dt + u*d(rho)/dt
% rhsE    = -rx.*(Dr*Ef)   + LIFT*(Fscale.*dEf/2); 
% % rhsp    = -rx.*(Dr*pf)   + LIFT*(Fscale.*dpf/2);
% 
% % [1       ]d(rho)/dt
% % [u  rho  ]d(u)dt
% % [       1]d(p)dt
% % Ainv = [1      0     0;
% %         -u/rho 1/rho 0]
% %         0      0     1]
% % rhsu = 
% 
% return

