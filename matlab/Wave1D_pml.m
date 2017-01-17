% wave

function [l2err]= Wave1D(smax,delta,p0sigma,p0approx)

if nargin < 2
    smax = 100;
    delta = .125;
    p0sigma = 0;
    p0approx = 0;
end

% Driver script for solving the 1D advection equations
Globals1D;

% Order of polymomials used for approximation
N = 4;

% Generate simple mesh
K1D = 32;
[Nv, VX, K, EToV] = MeshGen1D(-1.0,1.0,K1D);

% Initialize solver and construct grid and metric
StartUp1D;


% vmapP(1) = vmapM(end); % hack for periodic
% vmapP(end) = vmapM(1); % hack for periodic

% Set initial conditions
% p = cos(pi/2*x);
p = exp(-100*x.^2);
u = zeros(size(x));

global sigmax

M = 3;
L = 1-delta;
sigmax = (x > L).*(x-L).^M + (x < -L).*(-x-L).^M;

P1 = eye(N+1);

% sigmax = (x > L).*(smax./(delta - (x-L)) - smax/delta) + (x < -L).*(smax./(delta - (-x-L)) - smax/delta);
% sigmax(sigmax > 1e7) = 0;

if p0sigma
    for e = 1:K
        sigmax(:,e) = mean(sigmax(:,e));
    end
end

sigmax = smax*sigmax/max(sigmax(:));

if p0approx
    P1 = V*diag([1;zeros(N,1)])*inv(V); % reduce to const mode
end

sk = 1;
for e = 1:K
    if max(sigmax(:,e)) > 1e-8
        pmlK(sk) = e;
        sk = sk + 1;
    end
end
nonPML = setdiff(1:K,pmlK);

% sigmax(:) = smax;
% sigmax(:,nonPML) = 0;

% reverse profile
if 0
    sigmax = smax-sigmax; 
    sigmax(:,nonPML) = 0;
    title('reversed profile')
end

Mhat = inv(V*V');
Mp = Mhat*p*diag(J(1,:));
pK = p(:,nonPML); Mp = Mp(:,nonPML);
u0norm = sqrt(pK(:)'*Mp(:));

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
        [rhsp rhsu] = WaveRHS1D(p,u);
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        p = p + rk4b(INTRK)*resp;
        u = u + rk4b(INTRK)*resu;
        
        p(:,pmlK) = P1*p(:,pmlK);
        u(:,pmlK) = P1*u(:,pmlK);
    end;
    % Increment time
    time = time+dt;
    if nargin==0
        if 1 && mod(tstep,10)==0
            plot(x,p,'.-');
            axis([-1 1 -1 1])
            drawnow
        elseif 0% mod(tstep,10)==0
            Mp = Mhat*p*diag(J(1,:));
            pK = p(:,nonPML); Mp = Mp(:,nonPML);
            l2norm(sk) = sqrt(pK(:)'*Mp(:));
            tvec(sk) = time;
            semilogy(time,l2norm(sk),'.')
            hold on
            drawnow
            sk = sk + 1;
        end
    end
end;

Mp = Mhat*p*diag(J(1,:));
pK = p(:,nonPML); Mp = Mp(:,nonPML);
l2err = sqrt(pK(:)'*Mp(:)) / u0norm;

function [rhsp rhsu] = WaveRHS1D(p,u)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

global sigmax

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

% pml terms
rhsp = rhsp - sigmax.*p;
rhsu = rhsu - sigmax.*u;
% rhsp = -rhsp;
% rhsu = -rhsu;
return

