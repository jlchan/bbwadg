% wave

function Wave1D_SEM

% Driver script for solving the 1D advection equations
Globals1D;

% Order of polymomials used for approximation
N = 1;

CG = 1;

% Generate simple mesh
% K1D = 16;
K1D = 64;
[Nv, VX, K, EToV] = MeshGen1D(-1.0,1.0,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

% vmapP(1) = vmapM(end); % hack for periodic
% vmapP(end) = vmapM(1); % hack for periodic

rp = linspace(-1,1,25)';
Vp = Vandermonde1D(N,rp)/V;
xp = Vp*x;

% CG

pairs = sort([vmapM vmapP],2);
[~,id] = unique(pairs(:,1));
pairs = pairs(id,:);
gmap = zeros(Np,K);
for i = 1:size(pairs,1)
    gmap(pairs(i,:)) = i;
end
offset = size(pairs,1);
for e = 1:K
    gmap(2:N,e) = (1:N-1) + offset;
    offset = offset + N-1;
end
R = sparse(1:Np*K,gmap,1,Np*K,(K+1) + (N-1)*K);

% Set initial conditions
% p = cos(pi/2*x);
p = (x < 1/4).*(x > -1/4); 
pex = @(x) 1-(x > -1/3) + exp(-100*(x-1/3).^2);
p = pex(x);
% p = exp(-100*x.^2);
u = zeros(size(x));


% Solve Problem
FinalTime = 8;

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

M = kron(diag(J(1,:)),inv(V*V'));

sk = 1;
% outer time step loop
for tstep=1:Nsteps
    for INTRK = 1:5
        [rhsp rhsu] = WaveRHS1D(p,u);
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;               
        if CG
            resp(:) = R*((R'*M*R)\(R'*M*resp(:)));
            resu(:) = R*((R'*M*R)\(R'*M*resu(:)));
        end
        
        p = p + rk4b(INTRK)*resp;
        u = u + rk4b(INTRK)*resu;
        
        
    end;
    % Increment time
    time = time+dt;
    if 0
        if mod(tstep,10)==0
            plot(x,p,'o');
            hold on
            plot(xp,Vp*p,'-');
            axis([-1 1 -1.5 1.5])
            hold off
            drawnow       
        end
    end
end;

figure
% plot(x(:),p(:),'o','linewidth',2,'markersize',8,'MarkerFaceColor','g');
hold on
pp = Vp*p;
plot(xp(:),pp(:),'-.','linewidth',2);
plot(xp(:),pex(xp(:)),'-','linewidth',2);
axis([-1 1 -.5 1.5])
hold off
grid on
legend('FEM approximation','Exact solution')
set(gca,'fontsize',15)
print(gcf,'-dpng','CG.png')
% print(gcf,'-dpng','DG1.png')
% print(gcf,'-dpng','DG4.png')


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

return

