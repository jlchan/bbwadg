function [rho lam] = advec_spectra(NB,Ksub,K1D,smoothKnots,alphain)

Globals1D;

if nargin==0    
    NB = 6;
    Ksub = 32; %ceil(N/2);
    K1D = 1; 
    smoothKnots = 0;
end
global alpha
if nargin < 5
    alpha = 1;
else 
    alpha = alphain;
end
    
FinalTime = 1;

N = NB+Ksub-1; % number of sub-mesh splines

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

ndofs = K*(N+1); % total # dofs

% Initialize solver and construct grid and metric
StartUp1D;

vmapP(end) = vmapM(1);
vmapP(1) = vmapM(end); % make periodic
% assume uniform mesh for simplicity
rx = rx(1);

rp = linspace(-1,1,2000)';
Vp = Vandermonde1D(N,rp)/V;
xp = Vp*x;

%% make splines

rB = JacobiGL(0,0,N);
xB = Vandermonde1D(N,rB)/V * x;

[BVDM M Dr] = bsplineVDM(NB,Ksub,rB,smoothKnots); % VDM for interp, mass, M\S

[rq wq] = JacobiGQ(0,0,N);
[Bq] = bsplineVDM(NB,Ksub,rq,smoothKnots); % VDM for interp, mass, M\S
% M = Bq'*diag(wq)*Bq; % under-integrate mass matrix
[Vp] = bsplineVDM(NB,Ksub,rp,smoothKnots); % VDM for interp, mass, M\S

Vq = Vandermonde1D(N,rq)/V;
xq = Vq*x;

Mf = zeros(N+1,2);
Mf(1,1) = 1;
Mf(N+1,2) = 1;
LIFT = M\Mf;

% [Vp] = bsplineVDM(NB,Ksub,rp,smoothKnots);

global a
a = ones(N+1,K);
a(:,1) = 1;
% a(:,2) = 1;


%% check eigs

if 1
    A = zeros((N+1)*K);
    u = zeros(N+1,K);
    for i = 1:(N+1)*K
        u(i) = 1;
        rhsu = AdvecRHS1D(u);
%         keyboard
        A(:,i) = rhsu(:);
        u(i) = 0;
    end
    [W D] = eig(A);
    lam = diag(D);
 
    [rho id] = max(abs(lam));
    tau = alpha;
    hold on
    if abs(tau-1)<1e-8
        plot(lam,'o','linewidth',2,'markersize',8)
    elseif abs(tau)<1e-8
        plot(lam,'x','linewidth',2,'markersize',8)
    else
        plot(lam,'^','linewidth',2,'markersize',8)
    end
%     keyboard

plot(rp,Vp*W(:,id))

    return

%     keyboard
end

% h = legend('\tau = 0','\tau = 1/2','\tau = 1');set(h,'fontsize',14)
grid on
set(gca,'fontsize',14)
axis([-100 50 -110 110])

% % Set initial conditions
% % uex = @(x) exp(-10*sin(2*pi*x).^2);
% uex = @(x) sin(9*pi*x);
% uex = @(x) (x>-.75).*(x<-.25);
% % u = uex(x);
% % u = BVDM\uex(xB); % initialize w/lsq fit
% u = M\(Bq'*diag(wq)*uex(xq));
% 
% % plot(xp,uex(xp));return

function [rhsu] = AdvecRHS1D(u)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form field differences at faces
global a
global alpha
du = zeros(Nfp*Nfaces,K);
%du(:) = (a(vmapP).*u(vmapP)-a(vmapM).*u(vmapM)).*(nx(:)-(alpha)*abs(nx(:)))/2;
du(:) = (a(vmapP).*u(vmapP)-a(vmapM).*u(vmapM)).*nx(:)/2-(alpha)*1/2*(a(vmapP).*u(vmapP)-a(vmapM).*u(vmapM));
% du(:) = max(a(vmapP),a(vmapM)).*(u(vmapP)-u(vmapM)).*(nx(:)-(alpha)*abs(nx(:)))/2;
% du(:) = (u(vmapP)-u(vmapM)).*(nx(:)-(alpha)*abs(nx(:)))/2;

% compute right hand sides of the semi-discrete PDE
rhsu = -a.*rx.*(Dr*u) - LIFT*(Fscale.*(du));

return

