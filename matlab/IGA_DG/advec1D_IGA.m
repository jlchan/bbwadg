function [L2err] = advec1D_IGA(NB,Ksub,K1D,dt)

Globals1D;

if nargin==0    
    NB = 5;
    Ksub = 128; %ceil(N/2);
    K1D = 1; 
    if Ksub==1
        dt = 2/((NB+1)^2*K1D);
    else
        dt = 4/((NB+1)*Ksub*K1D);
    end
    smoothKnots = 0;
end

FinalTime = 1;

N = NB+Ksub-1; % number of sub-mesh splines

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

ndofs = K*(N+1) % total # dofs

% Initialize solver and construct grid and metric
StartUp1D;

vmapP(end) = vmapM(1);
vmapP(1) = vmapM(end); % make periodic
% assume uniform mesh for simplicity
rx = rx(1);

rp = linspace(-1,1,100)';
Vp = Vandermonde1D(N,rp)/V;
xp = Vp*x;

%% make splines

rB = JacobiGL(0,0,N);
xB = Vandermonde1D(N,rB)/V * x;

[BVDM M Dr] = bsplineVDM(NB,Ksub,rB,smoothKnots); % VDM for interp, mass, M\S

[rq wq] = JacobiGQ(0,0,N);
[Bq] = bsplineVDM(NB,Ksub,rq,smoothKnots); % VDM for interp, mass, M\S
% M = Bq'*diag(wq)*Bq; % under-integrate mass matrix

Vq = Vandermonde1D(N,rq)/V;
xq = Vq*x;

Mf = zeros(N+1,2);
Mf(1,1) = 1;
Mf(N+1,2) = 1;
LIFT = M\Mf;

[Vp] = bsplineVDM(NB,Ksub,rp,smoothKnots);

% Set initial conditions
% uex = @(x) exp(-10*sin(2*pi*x).^2);
uex = @(x) sin(9*pi*x);
uex = @(x) (x>-.75).*(x<-.25);
% u = uex(x);
% u = BVDM\uex(xB); % initialize w/lsq fit
u = M\(Bq'*diag(wq)*uex(xq));

% plot(xp,uex(xp));return
%% check eigs

if 1
    A = zeros((N+1)*K);
    u = zeros(N+1,K);
    for i = 1:(N+1)*K
        u(i) = 1;
        rhsu = AdvecRHS1D(u);
        A(:,i) = rhsu(:);
        u(i) = 0;
    end
    lam = eig(A);
    hold on
    plot(lam,'x')
    NB,max(abs(lam))
    max(real(lam))
    return

    keyboard
end

%%  Solve Problem

time = 0;

% Runge-Kutta residual storage
resu = zeros(N+1,K);

% xmin = min(abs(x(1,:)-x(2,:))); % compute time step size
% dt   = (.75/(K*N)); % should be O(h/N^2)
% dt   = .5*xmin; % should be O(h/N^2)

Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;

% outer time step loop
for tstep=1:Nsteps
    
    % low storage RK
    for INTRK = 1:5
        rhsu = AdvecRHS1D(u);
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u + rk4b(INTRK)*resu;        
    end;
    
    % Increment time
    time = time+dt;
    
    if nargin==0 && mod(tstep,10)==0 || tstep==Nsteps
        
        plot(xp,Vp*u,'-')
        hold on;
        plot(x,u,'o')
        plot(xp,uex(xp-time),'--');
        hold off
        axis([-1 1 -1 3])
        title(sprintf('Time = %f\n',time))
        drawnow
    end
end;

% err = uex(xp-time) - Vp*u;
% err = max(abs(err(:)))
[rq wq] = JacobiGQ(0,0,N+4);
Vq = Vandermonde1D(N,rq)/V;
xq = Vq*x;
[Bq] = bsplineVDM(NB,Ksub,rq,smoothKnots);
L2err = diag(wq)*abs(Bq*u-uex(xq-FinalTime)).^2*diag(J(1,:));
L2err = sqrt(sum(L2err(:)))


function [rhsu] = AdvecRHS1D(u)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form field differences at faces
alpha = 0;
du = zeros(Nfp*Nfaces,K);
du(:) = (u(vmapP)-u(vmapM)).*(nx(:)-(1-alpha)*abs(nx(:)))/2;

% compute right hand sides of the semi-discrete PDE
rhsu = -rx.*(Dr*u) - LIFT*(Fscale.*(du));

return

