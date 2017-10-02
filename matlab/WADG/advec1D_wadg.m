function advec1D_wadg

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

% vmapP(1) = vmapM(end); % make periodic
% vmapP(end) = vmapM(1);

rp = linspace(-1,1,50)';
Vp = Vandermonde1D(N,rp)/V;
xp=  Vp*x;

global uex Vq Pq

[rq wq] = JacobiGQ(0,0,2*N);
Vq = Vandermonde1D(N,rq)/V;
Pq = V*V'*Vq'*diag(wq);
xq =  Vq*x;

% Set initial conditions
d = 100;
uex = @(x) -1./(1 + exp(-d*(x-1/3))) + 1./(1 + exp(-d*(x+1/3))); % x > -1/3 & x < 1/3;%exp(-25*x.^2);
uex = @(x) (1+cos(pi*x)).*(abs(x)<1);
% uex = @(x) exp(-10.^2*(x+.125).^2);
% uex = @(x) -sin(pi*x);
% uex = @(x) 1 + (x < (-1+2/K1D));
% u = limit(uex(x));
u = uex(x);

% LIFT = diag(1./sum(inv(V*V'),2))*inv(V*V')*LIFT;

% Solve Problem
FinalTime = .25;

global a

a = @(x) 2+x;
% syms x(t) a
% eqn = diff(x,t,2) == 2+x;
% xSol(t) = dsolve(eqn)
% return

uexq = @(t) uex((2+xq)*exp(-t) - 2);
% uexq = uex(xq - FinalTime);
% plot(xq,uexq);return

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K);

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));

dt   = .25*xmin;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;


% outer time step loop
for tstep=1:Nsteps
    
    for INTRK = 1:5
        
        timeloc = time + rk4c(INTRK)*dt;
        [rhsu] = AdvecRHS1D(u);
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u + rk4b(INTRK)*resu;
    end
    
    %     % ssp-rk3
    %     ssprk = [1 .25 2/3];
    %     utmp = u;
    %     for i = 1:3
    %         [rhsu] = BurgersRHS1D(utmp);
    %         utmp = (1-ssprk(i))*u + ssprk(i)*(utmp + dt*rhsu);
    %     end
    %     u = utmp;
    
    
    if mod(tstep,5)==0
        
        plot(xp,Vp*u,'-','linewidth',2)
        hold on;
        plot(x,u,'o')
        %         plot(xe,VB\u,'o')
        
                plot(xq,uexq(time),'--','linewidth',2);
        hold off
        axis([-1 1 -1 3])
        title(sprintf('Time = %f\n',time))
        drawnow
    end
    % Increment time
    time = time+dt;
end;

plot(xq,uexq(FinalTime),'x','linewidth',2);
hold on;
plot(xp,Vp*u,'-','linewidth',2)
plot(x,u,'o')
% clf
% plot(xq,Vq*u-uexq,'x','linewidth',2);
% hold off
% axis([-1 1 -1 3])
title(sprintf('Time = %f\n',time))
legend('Exact sol at quad pts','Approximate solution')
set(gca,'fontsize',15)

% print(gcf,'-dpng','advecSol.png')
% keyboard


function [rhsu] = AdvecRHS1D(u)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;
global uex Pq Vq

% form field differences at faces
alpha = 0;
du = zeros(Nfp*Nfaces,K);

uM = u(vmapM);
uP = u(vmapP);
uP(1) = uex(-1);  
% uP(end) = uex(1); % BCs


global a
    aq = a(Vq*x);
    aM = a(x(vmapM));
    
if 0
    %     aq = 1;  aM = 1;

    du(:) = (uM-uP).*(aM.*nx(:)-(1-alpha)*abs(aM.*nx(:)))/2;
    rhsu = -Pq*((aq).*(Vq*(rx.*(Dr*u)))) + LIFT*(Fscale.*(du));
else
    du(:) = (uM-uP).*(nx(:)-(1-alpha)*abs(nx(:)))/2;
    rhsu = -rx.*(Dr*u) + LIFT*(Fscale.*(du));
    rhsu = Pq*(aq.*(Vq*rhsu));
end

return


