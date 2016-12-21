% wave

function [l2err]= Wave1D(smax,delta,p0sigma,p0approx)

% Driver script for solving the 1D advection equations
Globals1D;

% Order of polymomials used for approximation
N = 5;

% Generate simple mesh
a = .1;
K1D = 2*a/.0025;
[Nv, VX, K, EToV] = MeshGen1D(-a,a,K1D);

% Initialize solver and construct grid and metric
StartUp1D;


rp = linspace(-1,1,100);
Vp = Vandermonde1D(N,rp)/V;
xp = Vp*x;

% Set initial conditions
% p = cos(pi/2*x);
% p = exp(-100*x.^2);
p = zeros(size(x));
u = zeros(size(x));

% Solve Problem
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
% if (t < .5)
%     dp(1) = 2*(sin(pi*time/.5)-p(1));
%     du(1) = -2*u(1);
% end    

tau = 1;
pflux = .5*(tau*dp - du.*nx);
uflux = .5*(tau*du.*nx - dp).*nx;

% compute right hand sides of the semi-discrete PDE
f0 = 100; t0 = 1/f0;
rick = 1e4*(1- 2*(pi*f0*(t-t0))^2)*exp(-(pi*f0*(t-t0)).^2);
ptsrc = zeros(Np,K); 
ptsrc(round(Np*K/2)) = 1; 
ptsrc = V*V' * ptsrc;
rhsp = -rx.*(Dr*u) + LIFT*(Fscale.*pflux);
rhsu = -rx.*(Dr*p) + LIFT*(Fscale.*uflux) + rick*ptsrc;

return

