% wave

function [L2err dofs] = Wave1D_IGA(NBin,Ksubin,K1Din,smoothKnotsin)

% Driver script for solving the 1D advection equations
Globals1D;

if nargin==0
    NB = 3;
    Ksub = 16;
    K1D = 2;
    smoothKnots = 0;
else    
    NB = NBin;
    Ksub = Ksubin;
    K1D = K1Din;
    smoothKnots = smoothKnotsin;
end
N = NB+Ksub-1;

FinalTime = .5;

[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

rp = linspace(-1,1,100);
Vp = Vandermonde1D(N,rp)/V;
xp = Vp*x;

dofs = (N+1)*K1D;

%% make splines

rB = JacobiGL(0,0,N);
xB = Vandermonde1D(N,rB)/V * x;

[BVDM M Dr R rBq wBq Bq] = bsplineVDM(NB,Ksub,rB,smoothKnots); % VDM for interp, mass, M\S

% [rq wq] = JacobiGQ(0,0,N);
Vq = Vandermonde1D(N,rBq)/V;
xq = Vq*x; 
wJq = diag(wBq)*(Vq*J);

Vq = Bq;
% [Vq] = bsplineVDM(NB,Ksub,rq,smoothKnots); % VDM for interp, mass, M\S
% M = Bq'*diag(wq)*Bq; % under-integrate mass matrix

Mf = zeros(N+1,2);
Mf(1,1) = 1;
Mf(N+1,2) = 1;
% keyboard
LIFT = M\Mf;

[Vp] = bsplineVDM(NB,Ksub,rp,smoothKnots);

Pq = M\(Vq'*diag(wBq));

%% Set initial conditions

k = 3;
pex = @(x,t) cos(k*pi/2*x).*cos(k*pi*t/2);
p = Pq*pex(xq,0); 
% p = exp(-100*x.^2);
u = zeros(size(x));

%% compute eigs
if 1
    U = zeros(Np*K,2);
    A = zeros(2*Np*K);
    for i = 1:2*Np*K
        U(i) = 1;
        [rhsp rhsu] = WaveRHS1D(reshape(U(:,1),Np,K),reshape(U(:,2),Np,K),0);
        A(:,i) = [rhsp(:);rhsu(:)];
        U(i) = 0;
    end
    
    lam = eig(A);
    rho = max(abs(lam))
    return
end
%% Solve Problem

time = 0;

% Runge-Kutta residual storage
resp = zeros(Np,K);
resu = zeros(Np,K);

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
dt  = .125*xmin;
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
    
    if nargin==0 && mod(tstep,10)==0
        clf
        pp = Vp*p;
        %semilogy(xp,abs(pex(xp,time)-pp),'-');
         plot(xp,pp,'r-');
%         hold on;
%         up = Vp*u;
%         plot(xp,up,'r-');
        axis([-1 1 -1 1])
        drawnow
    end
end;

L2err = sqrt(sum(sum(wJq.*(Vq*p - pex(xq,FinalTime)).^2)));

function [rhsp rhsu] = WaveRHS1D(p,u,t)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form field differences at faces
dp = zeros(Nfp*Nfaces,K);
du = zeros(Nfp*Nfaces,K);
dp(:) = p(vmapP)-p(vmapM);
du(:) = u(vmapP)-u(vmapM);

dp(mapB) = -2*p(vmapB);
% du(mapB) = -2*u(vmapB);

tau = 1;
pflux = .5*(tau*dp - du.*nx);
uflux = .5*(tau*du.*nx - dp).*nx;

% compute right hand sides of the semi-discrete PDE
rhsp = -rx.*(Dr*u) + LIFT*(Fscale.*pflux);
rhsu = -rx.*(Dr*p) + LIFT*(Fscale.*uflux);

return

