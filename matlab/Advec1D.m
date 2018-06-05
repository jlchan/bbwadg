% wave

% function Advec1D

% Driver script for solving the 1D advection equations
Globals1D;

% Order of polymomials used for approximation
N = 3;
K1D = 32;
CFL = .5;
FinalTime = 10;

[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
StartUp1D;

vmapP(1) = vmapM(end); % hack for periodic
vmapP(end) = vmapM(1); % hack for periodic

global a
a = 1;
global tau
tau = 1;

rp = linspace(-1,1,100);
Vp = Vandermonde1D(N,rp)/V;
xp = Vp*x;

% Set initial conditions
uex = @(x) exp(-100*x.^2);
% uex = @(x) sin(4*pi*x);
% uex = @(x) (abs(x) < 1/2);
[rq wq] = JacobiGQ(0,0,N+1);
Vq = Vandermonde1D(N,rq)/V;
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));
xq = Vq*x;
wJq = diag(wq)*(Vq*J);

u = Pq*uex(xq);


% Solve Problem
time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K);

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
dt  = CFL*xmin;
% dt = min(dt,1/smax);
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

sk = 1;
% outer time step loop
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        [rhsu] = AdvecRHS1D(u,timelocal);
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u + rk4b(INTRK)*resu;        
    end;
    % Increment time
    time = time+dt;
    
    if 1 && mod(tstep,50)==0 || tstep==Nsteps
        clf
        pp = Vp*u;
        plot(xp,pp,'r--','linewidth',2);
        hold on
        title(sprintf('time = %f\n',time));
        %         plot(x,u,'bo','linewidth',2,'markersize',8,'MarkerFaceColor',[.49 1 .63]);
        axis([-1,1,-2 2])
        drawnow
    end
    
    energy(tstep) = sum(sum(wJq.*(Vq*u)));
end;

clf
plot(xp,uex(xp),'b-','linewidth',2);
hold on
pp = Vp*u;
plot(x([1 N+1],:),u([1 N+1],:),'ro','linewidth',2,'markersize',8,'MarkerFaceColor',[.49 1 .63])
plot(xp,pp,'r--','linewidth',2);
% axis([-1,1,-.25,1.25])
set(gca,'fontsize',15)
grid on
% 
% figure(2)
% semilogy((1:Nsteps)*dt,energy,'--','linewidth',2)
% hold on
% set(gca,'fontsize',15)
% grid on
% xlabel('Time','fontsize',15)
% ylabel('$\|u\|^2$','fontsize',18,'Interpreter','LaTeX')



function [rhsu] = AdvecRHS1D(u,time)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form field differences at faces
global a
global tau
du = zeros(Nfp*Nfaces,K);
uP = reshape(u(vmapP),2,K);
% uP(1) = exp(-100*(time-.25).^2);
% du(:) = (uP-u(vmapM)).*(a*nx(:)-(tau)*abs(a*nx(:)))/2;
uM = reshape(u(vmapM),2,K);
du = a*.5*(uP-uM).*nx - .5*abs(a)*tau*(uP-uM);

% compute right hand sides of the semi-discrete PDE
rhsu = -a*rx.*(Dr*u) - LIFT*(Fscale.*(du));
end


