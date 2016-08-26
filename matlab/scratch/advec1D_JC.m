function advec1D_JC

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;

% Order of polymomials used for approximation
N = 3;

% Generate simple mesh
K1D = 32;
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

vmapP(1) = vmapM(end); % make periodic
% vmapP(1) = vmapM(end); % hack for periodic
% vmapP(end) = vmapM(1); % hack for periodic

% Set initial conditions
d = 100;
u = -1./(1 + exp(-d*(x-1/3))) + 1./(1 + exp(-d*(x+1/3)));%x > -1/3 & x < 1/3;%exp(-25*x.^2);
% u = -sin(pi*x);

% Solve Problem
FinalTime = 2;

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K);

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));

dt   = .25*xmin;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;

% advection speed
a = 1;

% outer time step loop
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        [rhsu] = AdvecRHS1D(u, timelocal, a);
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u+rk4b(INTRK)*resu;
    end;
    if mod(tstep,10)
        plot(x,u,'o-')
        axis([-1 1 -2 2])
        drawnow
    end
    % Increment time
    time = time+dt;
end;


function [rhsu] = AdvecRHS1D(u,time, a)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form field differences at faces
alpha = 1;
du = zeros(Nfp*Nfaces,K);
du(:) = (u(vmapM)-u(vmapP)).*(a*nx(:)-(1-alpha)*abs(a*nx(:)))/2;
% a = u; du(:) = (u(vmapM)-u(vmapP)).*(a(vmapM).*nx(:)-(1-alpha)*abs(a(vmapM).*nx(:)))/2;

% % impose boundary condition at x=0
% uin = -sin(a*time);
% du (mapI) = (u(vmapI)- uin ).*(a*nx(mapI)-(1-alpha)*abs(a*nx(mapI)))/2;
% du (mapO) = 0;

% compute right hand sides of the semi-discrete PDE
rhsu = -a.*rx.*(Dr*u) + LIFT*(Fscale.*(du));
return
