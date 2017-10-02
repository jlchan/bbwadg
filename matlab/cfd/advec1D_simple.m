function advec1D_simple

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u
Globals1D;

% Order of polymomials used for approximation
N = 4;
Ksub = 8; %ceil(N/2);
NB = N+Ksub-1; % number of sub-mesh splines

% Generate simple mesh
K1D = 4; 
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

K*(NB+1) % total # dofs

% Initialize solver and construct grid and metric
StartUp1D;

rp = linspace(-1,1,50)';
Vp = Vandermonde1D(N,rp)/V;
xp = Vp*x;

% make splines
rB = JacobiGL(0,0,NB);
VB = Vandermonde1D(NB,rB);
xB = Vandermonde1D(N,rB)/V * x;

[BVDM M Dr] = bsplineVDM(N,Ksub,rB);
Mf = zeros(NB+1,2);
Mf(1,1) = 1;
Mf(NB+1,2) = 1;
LIFT = M\Mf;

[Vp] = bsplineVDM(N,Ksub,rp);

% plot(rp,Vp,'--'); return
vmapM = reshape(vmapM,2,K);
vmapP = reshape(vmapP,2,K);
for e = 1:K
    vmapM(:,e) = [1; NB+1] + (e-1)*(NB+1);
    vmapP(:,e) = [0; NB+2] + (e-1)*(NB+1);
end
vmapM = vmapM(:);
vmapP = vmapP(:);

vmapP(end) = vmapM(1);
vmapP(1) = vmapM(end); % make periodic

% assume uniform mesh for simplicity
rx = rx(1);

% Set initial conditions
% uex = @(x) exp(-25*sin(2*pi*x).^2);
uex = @(x) sin(2*pi*x);
% u = uex(x);
u = BVDM\uex(xB); % initialize w/lsq fit

%% check eigs

if 1
    A = zeros((NB+1)*K);
    u = zeros(NB+1,K);
    for i = 1:(NB+1)*K
        u(i) = 1;
        rhsu = AdvecRHS1D(u);
        A(:,i) = rhsu(:);
        u(i) = 0;
    end
    lam = eig(A);
    hold on
    plot(lam,'x')
    NB,max(abs(lam))
    return
    % keyboard
end

%%

% Solve Problem
FinalTime = 1;

time = 0;

% Runge-Kutta residual storage
resu = zeros(NB+1,K);
xmin = min(abs(x(1,:)-x(2,:))); % compute time step size
dt   = (.75/(K*NB)); % should be O(h/N^2)
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
    
    if mod(tstep,10)==0
        
        plot(xp,Vp*u,'-')
        hold on;
%         plot(x,u,'o')
        plot(xp,uex(xp-time),'--');
        hold off
        axis([-1 1 -1 3])
        title(sprintf('Time = %f\n',time))
        drawnow
    end
end;

% err = uex(xp-time) - Vp*u;
% err = max(abs(err(:)))
[rq wq] = JacobiGQ(0,0,NB+4);
Vq = Vandermonde1D(N,rq)/V;
xq = Vq*x;
[Bq] = bsplineVDM(N,Ksub,rq);
L2err = diag(wq)*abs(Bq*u-uex(xq-FinalTime)).^2*diag(J(1,:));
L2err = sqrt(sum(L2err(:)))


function [rhsu] = AdvecRHS1D(u)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form field differences at faces
alpha = 0;
du = zeros(Nfp*Nfaces,K);
du(:) = (u(vmapM)-u(vmapP)).*(nx(:)-(1-alpha)*abs(nx(:)))/2;

% compute right hand sides of the semi-discrete PDE
rhsu = -rx.*(Dr*u) + LIFT*(Fscale.*(du));
return

