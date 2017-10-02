function advec1D_bern_artvisc

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;

% Order of polymomials used for approximation
N = 7;

% Generate simple mesh
K1D = 16;
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

vmapP(1) = vmapM(end); % make periodic
vmapP(end) = vmapM(1);

rp = linspace(-1,1,25)';
Vp = Vandermonde1D(N,rp)/V;
xp=  Vp*x;

global VB DB VBe W A
W = get_BB_smoother(N);

re = linspace(-1,1,N+1)';
[VB VBr] = bern_basis_1D(N,r);
DB = VB\VBr;
DB(abs(DB)<1e-8) = 0;
VBe = bern_basis_1D(N,re);
xe = Vandermonde1D(N,re)/V * x;

% Set initial conditions
d = 100;
uex = @(x) -1./(1 + exp(-d*(x-1/3))) + 1./(1 + exp(-d*(x+1/3)));%x > -1/3 & x < 1/3;%exp(-25*x.^2);
uex = @(x) (x > -3/4 & x < -1/4) + exp(-50*(x-1/2).^2);
% uex = @(x) 1-cos(pi*x);
% uex = @(x) sin(pi*x);
% u = limit(uex(x));
% u = uex(xe);
u = uex(x);

A = diag(ones(N+1,1)) - diag(ones(N,1),1);

% Solve Problem
FinalTime = 2;

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K);

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));

dt   = .4*xmin;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;
ssprk = [1 .25 2/3]; % ssp coeffs

% outer time step loop
for tstep=1:Nsteps
    
    % ssp-rk3      
    utmp = u;
    for i = 1:3        
        [TV Klim] = limit(utmp);        
        [rhsu] = AdvecRHS1D(utmp,A,Klim);
        plot(xp,Vp*u,'-')
        utmp = (1-ssprk(i))*u + ssprk(i)*(utmp + dt*rhsu);        
    end    
    u = utmp;

    % Increment time
    time = time+dt;
    
    if mod(tstep,5)==0
        
        plot(xp,Vp*u,'-')
        %         plot(x,VB*uex(xe-time),'.-')
        hold on;
        plot(x,u,'o')
        plot(xe,VB\u,'o')
        
        plot(xp,uex(xp-time),'--');
        plot(x(:,Klim),ones(size(x(:,Klim))),'x')
        hold off
        axis([-1 1 -1 3])
        title(sprintf('Time = %f\n',time))
        drawnow
    end
    
end;

plot(xp,Vp*u,'-')
hold on;
plot(x,u,'o')
plot(xp,uex(xp-time),'--');
axis([-1 1 -1 3])
% keyboard

function [rhsu] = AdvecRHS1D(u,A,Klim)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

global VB

% form field differences at faces
alpha = 0;
du = zeros(Nfp*Nfaces,K);
du(:) = (u(vmapM)-u(vmapP)).*(nx(:)-(1-alpha)*abs(nx(:)))/2;

% compute right hand sides of the semi-discrete PDE
rhsu = -rx.*(Dr*u) + LIFT*(Fscale.*(du));

uB = VB\u(:,Klim);
TV = A*uB;
ep = sum(abs(TV));

% % if length(ep)>0
% %     ep = ep/max(abs(ep(:)));
% % end
% TV = TV.*abs(TV);
% uav = A'*TV;
% % uav = A'*A*uB;
% uav = uav*diag(ep);

% uav = A'*A*uB*diag(ep);
% uav = VB*(uav);
% rhsu(:,Klim) = rhsu(:,Klim) - uav;
return


function [TV Klim] = limit(u)

Globals1D
global VB DB VBe W

% convert to bernstein
uB = VB\u;

% find K with large TV
TV = 0;
for i = 1:N  
    TV = TV + abs(uB(i,:) - uB(i+1,:));
end
% TV = TV./(N*max(abs(uB(:))));
Klim = find(TV > (N+1));
TV = TV(:,Klim);
