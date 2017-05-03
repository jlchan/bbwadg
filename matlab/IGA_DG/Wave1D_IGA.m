% wave

function rho = Wave1D_IGA(NB,Ksub,K1D)

% Driver script for solving the 1D advection equations
Globals1D;

if nargin==0
    N = 5;
    K1D = 32;
else
    
end

[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

rp = linspace(-1,1,100);
Vp = Vandermonde1D(N,rp)/V;
xp = Vp*x;

% Set initial conditions
% p = cos(pi/2*x);
p = exp(-100*x.^2);
u = zeros(size(x));

%% compute eigs

U = zeros(Np*K,2);
A = zeros(2*Np*K);
for i = 1:2*Np*K
    U(i) = 1;
    [rhsp rhsu] = WaveRHS1D(reshape(U(:,1),Np,K),reshape(U(:,2),Np,K),0);
    A(:,i) = [rhsp(:);rhsu(:)];
    U(i) = 0;   
end

lam = eig(A);
rho = max(abs(lam));
return
%% Solve Problem
FinalTime = 2;

time = 0;

% Runge-Kutta residual storage
resp = zeros(Np,K);
resu = zeros(Np,K);

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
dt  = .5*xmin;
% dt = min(dt,1/smax);
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

sk = 1;
% outer time step loop
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        [rhsp rhsu] = WaveRHS1D(p,u,timelocal);
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        p = p + rk4b(INTRK)*resp;
        u = u + rk4b(INTRK)*resu;
    end;
    % Increment time
    time = time+dt;
    
    if 1 && mod(tstep,10)==0
        clf
        pp = Vp*p;
        plot(xp,pp,'-');
%         hold on;
%         up = Vp*u;
%         plot(xp,up,'r-');
%         axis([-1 1 -1 1])
        drawnow
    end
end;


function [rhsp rhsu] = WaveRHS1D(p,u,t)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form field differences at faces
dp = zeros(Nfp*Nfaces,K);
du = zeros(Nfp*Nfaces,K);
dp(:) = p(vmapP)-p(vmapM);
du(:) = u(vmapP)-u(vmapM);

dp(mapB) = 0;
du(mapB) = -2*u(vmapB);

tau = 1;
pflux = .5*(tau*dp - du.*nx);
uflux = .5*(tau*du.*nx - dp).*nx;

% compute right hand sides of the semi-discrete PDE
rhsp = -rx.*(Dr*u) + LIFT*(Fscale.*pflux);
rhsu = -rx.*(Dr*p) + LIFT*(Fscale.*uflux);

return

