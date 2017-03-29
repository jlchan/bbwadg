function Euler1D

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;

% Order of polymomials used for approximation
N = 4;
K1D = 16;
FinalTime = 5;

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

W = get_BB_smoother(N);
W = VBe;
W = VB*W/VB;

% Set initial conditions
global gamma
gamma = 1.4;

[rho m E] = vortexSolution(x,0);

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

dt   = .2*xmin;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;

M = inv(V*V');

% outer time step loop
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
        
        [rho m E] = limit(rho,m,E);        
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
        
        ylim([-1 8])
        title(sprintf('Time = %f\n',time))
        drawnow
    end
    
    enorm(tstep) = sqrt(sum(sum(J.*(M*(rho.^2 + m.^2 + E.^2)))));
    
    if mod(tstep,1000)==0
        disp(sprintf('tstep = %d out of %d\n',tstep,Nsteps))
    end
end;
figure;
plot(enorm)

function [rho m E] = limit(rho,m,E)

Globals1D
global VB W
TV = TV1D(N,VB\rho) + TV1D(N,VB\m) + TV1D(N,VB\E);
ids = find(TV > 2*(N+1));
rho(:,ids) = W*rho(:,ids);
m(:,ids)   = W*m(:,ids);
E(:,ids)   = W*E(:,ids);

function TV = TV1D(N,u)

TV = 0;
for i = 1:N
    TV = TV + abs(u(i,:) - u(i+1,:));
end

function [rho m E] = vortexSolution(x,t)

global gamma

x0 = 5; 
% r0 = .25;
% rho0 = 1;
% p0 = 1/gamma; 
% c0 = sqrt(gamma*p0./rho0);
% r2 = (x-x0-t).^2;
% PP = exp(.5*(1-r2/r0));
% rho = 1 - .5*(gamma-1)*PP.^2;
% rho = rho.^(1/(gamma-1));
% p = rho.^gamma;
% u = 1-c0*PP;
% m = rho.*u;
% E = p/(gamma-1) + .5*rho.*u.^2;

% rho = 1+exp(-(x-x0).^2); 
% % rho = 1 + (x > x0);
% p = rho.^gamma;
% u = 0*rho;
% E = p/(gamma-1) + .5*rho.*u.^2;
% m = rho.*u;

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

% compute maximum velocity for LF flux
u = m./rho;
pres = (gamma-1.0)*(E - 0.5*rho.*u.^2);
cvel = sqrt(gamma*pres./rho); 
lm   = abs(u) + cvel;

% Compute fluxes
rhof = m; 
mf   = rho.*u.^2 + pres; 
Ef   = (E+pres).*u;

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


