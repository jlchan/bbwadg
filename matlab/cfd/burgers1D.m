function L2err = burgers1D(Nin,K1D)
% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;

if nargin==0
    N = 7;
    K1D = 16;
else
    N = Nin;
end
CFL = .125;
FT = 1; % 1.06765;
global tau
tau = 1;
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% VX(2:end-1) = VX(2:end-1) + 1/K1D*randn;

% Initialize solver and construct grid and metric
StartUp1D;

% make periodic
vmapP(1) = vmapM(end); 
vmapP(end) = vmapM(1);

rp = linspace(-1+1e-4,1-1e-4,50)';
Vp = Vandermonde1D(N,rp)/V;
xp=  Vp*x;

global uex Vq Pq Prq shockSol Vf

Vf = Vandermonde1D(N,[-1;1])/V;
 
% Nq = ceil((3*N-1)/2); % exact integration
% Nq = 2*N+2; % same as overkill integration
Nq = N; % inexact integration
[rq wq] = JacobiGL(0,0,Nq);

Vq = Vandermonde1D(N,rq)/V;
Vrq = GradVandermonde1D(N,rq)/V;
Pq = V*V'*Vq'*diag(wq);
Prq = V*V'*Vrq'*diag(wq);
xq =  Vq*x;

global M
M = inv(V*V');
% LIFT = diag(1./sum(M,2))*(M*LIFT);
% keyboard

shockSol = 1;
if shockSol
    uex = @(x) 1 + (x < (-1 + .25));
    uex = @(x) -1 + (x>0)*2;    
    FinalTime = 1/2;
else
    uex = @(x) -sin(pi*x);
    FinalTime = FT;
    
end
u = Pq*uex(xq);

M = inv(V*V');
invM = V*V';


% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
dt   = CFL*xmin;
dt   =  0.0036;
resu = zeros(Np,K);
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;

D = diag([ones(N,1); 1]);
F = V*D/V;

% outer time step loop
tvec = round(Nsteps/FinalTime * .25 *[1:7]);
figure(1)
for tstep=1:Nsteps
    
    % low storage RK
    for INTRK = 1:5
        rhsu = BurgersRHS1D(u);
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u + rk4b(INTRK)*resu;
        
%         u = F*u;
%         u = limit1D(u);
    end;
    
    time = (tstep-1)*dt;
    
    if nargin==0 && (mod(tstep,5)==0 || tstep==Nsteps)
        
        plot(xp,Vp*u,'b-','linewidth',2)
        if shockSol
            hold on;
            plot(xp,uex(xp-1.5*time));
            hold off
        end
        axis([-1 1 -4 4])
        set(gca,'fontsize',15);grid on
        title(sprintf('Time = %f\n',tstep*dt),'fontsize',15)
        
        drawnow
        
[val id] = min(abs(tstep-tvec));
        if val<1e-8
%             tag = '~/Desktop/bbwadg/talks/TAMES2017/figs/burgersStableEC_%d.png';
%             tag = '~/Desktop/bbwadg/talks/TAMES2017/figs/burgersStableLF_%d.png';
%             tag = '~/Desktop/bbwadg/talks/TAMES2017/figs/burgersStable_%d.png';
%             tag = '~/Desktop/bbwadg/talks/TAMES2017/figs/burgersUnstable_%d.png';
%             print(gcf,'-dpng',sprintf(tag,id))            
%             keyboard
%             pause(.1)
        end
    end
    
    energy(tstep) = sum(sum(diag(wq)*(Vq*u).^2*J(1)));
    
end;
return
% print(gcf,'-dpng',sprintf(tag,length(tvec)+1))
% print(gcf,'-dpng',sprintf('~/Desktop/bbwadg/talks/TAMES2017/figs/burgersStable_%d.png',4))
% print(gcf,'-dpng',sprintf('~/Desktop/bbwadg/talks/TAMES2017/figs/burgersUnstable_%d.png',4))
figure(2)
plot(dt*(1:Nsteps),energy,'-','linewidth',2)
set(gca,'fontsize',15)
grid on
xlabel('Time','fontsize',15)
ylabel('$\|u\|^2$','fontsize',18,'Interpreter','LaTeX')
%ytickformat('%5.5f')
ytickformat('%2.1f')
ylim([0 1.1])
print(gcf,'-dpng','~/Desktop/bbwadg/talks/TAMES2017/figs/burgersSplitEnergyEC.png')
% print(gcf,'-dpng','~/Desktop/bbwadg/talks/TAMES2017/figs/burgersSplitEnergyLF.png')
        



function [rhsu] = BurgersRHS1D(u)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;
global uex Pq Vq Prq shockSol Vf
global tau
% form field differences at faces
alpha = 0;
du = zeros(Nfp*Nfaces,K);
uavg = zeros(Nfp*Nfaces,K);
uM = reshape(u(vmapM),Nfp*Nfaces,K);
uP = reshape(u(vmapP),Nfp*Nfaces,K);

if shockSol
    uP(1) = uex(-1);
    uP(end) = uex(1); % BCs
    f1 = uex(-1).^2;
    fend = uex(1).^2;
end

opt = 3;
if opt==1 % conservative
    du(:) = .5*(uP-uM).*(uM.*nx - abs(uM.*nx));
    
    rhsu = -Pq*((Vq*u).*(Vq*(rx.*(Dr*u)))) - LIFT*(Fscale.*(du)); % non-conservative form
%     rhsu = -rx.*(Dr*(Pq*.5*(Vq*u).^2)) - LIFT*(Fscale.*du); % quadrature projection (non-conservative!) form
%     rhsu = -rx.*(Dr*(.5*u.^2)) - LIFT*(Fscale.*du); % quadrature conservation form
    
%     % entropy-conserving flux (strong form)
%     fq = .5*(Vq*u).^2;        
%     f = Pq*fq; 
%     fM = reshape(f(vmapM),Nfp*Nfaces,K);   
%     du(:) = -(uP.^2 + uP.*uM + uM.^2)/6.*nx + fM.*nx + .5*abs(uM).*(uP-uM);
%     rhsu = -rx.*(Dr*f) + LIFT*(Fscale.*du);
            

elseif opt==2 % discrete split
    
    du(:) = uP-uM;
    Du = rx.*(Dr*u) + .5*LIFT*(Fscale.*du.*nx); % Dh*f
    uDu = Pq*((Vq*u).*(Vq*Du));
    
    f = Pq*((Vq*u).^2); % project flux
    fP = f(vmapP); 
    fM = f(vmapM);
    du(:) = fP - fM;
    
    if shockSol
        % shock BCs
        du(1) = f1 - fM(1);
        du(end) = fend - fM(end);
    end
    Df = rx.*(Dr*f) + .5*LIFT*(Fscale.*du.*nx); % Dh*f
    rhsu = -(uDu + Df)*1/3;
    
    % add lax friedrichs
    du = tau*.5*(uP-uM).*abs(uM);
    rhsu = rhsu + LIFT*(Fscale.*du);

elseif opt==3 % continuous split + discretization
    
    du(:) = uM.*(uP-uM); 
    uDu = Pq*((Vq*u).*(Vq*(rx.*(Dr*u)))) + .5*LIFT*(Fscale.*du.*nx); % Dh*f    
        
    uavg(:) = .5*(uP.^2 + uM.^2);
    %Df = -rx.*(Prq*(Vq*u).^2) + LIFT*(Fscale.*uavg.*nx); % Dh*f
    Df = rx.*(Dr*Pq*(Vq*u).^2) + LIFT*(Fscale.*(uavg - Vf*(Pq*(Vq*u).^2)).*nx); % Dh*f
    
    rhsu = -(uDu + Df)*1/3;    
    
    % lax friedrichs
    du = tau*(uP-uM).*abs(uM);
    rhsu = rhsu + .5*LIFT*(Fscale.*du);
    
elseif opt==4 % discrete split + interpolation instead of quadrature    
    
    du(:) = uP-uM;
    Du = rx.*(Dr*u) + .5*LIFT*(Fscale.*du.*nx); % Dh*f
    global M
    uDu = M\(u.*(M*Du)); % energy stable, may not be accurate in curvilinear coordinates though
    
    f = u.^2; % project interpolated flux
    fP = f(vmapP); 
    fM = f(vmapM);
    du(:) = fP - fM;
    
    if shockSol
        % shock BCs
        du(1) = f1 - fM(1);
        du(end) = fend - fM(end);
    end
    Df = rx.*(Dr*f) + .5*LIFT*(Fscale.*du.*nx); % Dh*f
    rhsu = -(uDu + Df)*1/3;
    
    % add lax friedrichs
    du = tau*.5*(uP-uM).*abs(uM);
    rhsu = rhsu + LIFT*(Fscale.*du);

elseif opt==5 % non-conservative part
    
    du(:) = uM.*(uP-uM); 
    uDu = Pq*((Vq*u).*(Vq*(rx.*(Dr*u)))) + .5*LIFT*(Fscale.*du.*nx); % Dh*f    
            
    rhsu = -uDu;    
    
    % lax friedrichs
    du = (uP-uM).*abs(uM);
    rhsu = rhsu + .5*LIFT*(Fscale.*du);
    
end

return


