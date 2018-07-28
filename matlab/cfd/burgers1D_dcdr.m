Globals1D;

N = 4;
K1D = 16;
CFL = .125;
FinalTime = 1;
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

global uex Vq Pq Prq shockSol Vf DNr Lq

Nq = N; % inexact integration
[rq wq] = JacobiGQ(0,0,Nq);

Vf = Vandermonde1D(N,[-1;1])/V;

% Nq = ceil((3*N-1)/2); % exact integration
% Nq = 2*N+2; % same as overkill integration

Vq = Vandermonde1D(N,rq)/V;
Vrq = GradVandermonde1D(N,rq)/V;
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq));
xq =  Vq*x;

% LIFT = diag(1./sum(M,2))*(M*LIFT);
wf = [1;1];
Lq = M\(Vf'*diag(wf));

WN = diag([wq;wf]);
DNr = [Vq*Dr*Pq - .5*Vq*Lq*diag([-1,1])*Vf*Pq .5*Vq*Lq*diag([-1,1]);
    -.5*diag([-1,1])*Vf*Pq 0*.5*diag([-1,1])];

%%


u = Pq*uex(xq);
% plot(x,u,'o')
% return

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
dt   = CFL*xmin;
resu = zeros(Np,K);
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;


% outer time step loop
figure(1)
for tstep=1:Nsteps
    
    % low storage RK
    for INTRK = 1:5
        rhsu = BurgersRHS1D(u);
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u + rk4b(INTRK)*resu;

    end;
    
    time = (tstep-1)*dt;
    
    if  (mod(tstep,5)==0 || tstep==Nsteps)
        
        plot(xp,Vp*u,'b-','linewidth',2)        
        axis([-1 1 -4 4])
        set(gca,'fontsize',15);grid on
        title(sprintf('Time = %f\n',tstep*dt),'fontsize',15)
        
        drawnow                
    end
    
    energy(tstep) = sum(sum(diag(wq)*(Vq*u).^2*J(1)));
    
end;

return



function [rhsu] = BurgersRHS1D(u)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;
global tau DNr Vq Vf Lq Pq

% form field differences at faces
uM = reshape(u(vmapM),Nfp*Nfaces,K);
uP = reshape(u(vmapP),Nfp*Nfaces,K);

uN = [Vq;Vf]*u;
fS = @(uL,uR) (1/6)*(uL.^2 + uL.*uR + uR.^2);

fSf = fS(uM,uP);
rhsu = zeros(N+1,K);
for e = 1:K    
    [ui, uj] = meshgrid(uN(:,e));    
    rhsu(:,e) = [Pq Lq]*sum(DNr.*fS(ui,uj),2);
end

rhsu = -(2*(rx.*rhsu) + Lq*(Fscale.*nx.*fSf));

end


