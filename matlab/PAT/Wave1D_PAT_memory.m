% wave

function Wave1D_simple

% Driver script for solving the 1D advection equations
Globals1D;

% Order of polymomials used for approximation
N = 4;
FinalTime = 8;
K1D = 8;

global tau0; tau0 = 1;

plotSol = 1;
useTraces = 1;

maxiter = 50;

bcOpt = 'rc'; %rd for dirichlet, rc for characteristic vars

[Nv, VX, K, EToV] = MeshGen1D(-1.0,1.0,K1D);


% Initialize solver and construct grid and metric
StartUp1D;

rp = linspace(-1,1,25)';
Vp = Vandermonde1D(N,rp)/V;
xp = Vp*x;

global Vq Pq c2 
[rq wq] = JacobiGQ(0,0,N+1);
Vq = Vandermonde1D(N,rq)/V;
xq = Vq*x;
Pq = (V*V') * (Vq'*diag(wq));
c2 = 1 + 0*exp(-25*xq.^2);
c2(:,1) = 1; % make last elem homogeneous for ABCs
c2(:,end) = 1; % make last elem homogeneous for ABCs

global kappa0
kappa0 = 0*exp(-25*xq.^2);

% Set initial conditions
% p = cos(pi/2*x);
% pex = @(x) 1-(x > -1/3) + exp(-100*(x-1/3).^2);
pex = @(x) (abs(x + 1/2)<.25) + exp(-5^2*(x-1/2).^2);
% plot(xp,pex(xp));return
pex = @(x) exp(-16*(x-.125).^2);

% pex = @(x) (abs(x)<.5);

p = Pq*pex(xq);
u = zeros(size(x));

p_init = p;

%% memory term

qfun = @(x) ones(size(x));
Phi = zeros(size(x)); % initalize to average of p assuming q(x) = 1

wJq = diag(wq)*(Vq*J);
%% time stepping

% Runge-Kutta residual storage
resp = zeros(Np,K);
resu = zeros(Np,K);
resphi = zeros(Np,K);

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
dt  = .5*xmin;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

M = kron(diag(J(1,:)),inv(V*V'));

% outer time step loop
pb = zeros(length(mapB),Nsteps*5);
ub = zeros(length(mapB),Nsteps*5);
for tstep=1:Nsteps
    for INTRK = 1:5
        %[rhsp rhsu rhsPhi] = WaveRHS1D(p,u,Phi,dt,'abc',[],[]);
        [rhsp rhsu rhsPhi] = WaveRHS1D(p,u,Phi,dt,'d',[],[]);
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;               
        resphi = rk4a(INTRK)*resphi + dt*rhsPhi;
        
        p = p + rk4b(INTRK)*resp;
        u = u + rk4b(INTRK)*resu;                
        Phi = Phi + rk4b(INTRK)*resphi;
        
        id = Nsteps*5 - ((tstep-1)*5 + INTRK) + 1;        
        pb(:,id) = p(vmapB);        
        ub(:,id) = u(vmapB);        
    end;
        
    if plotSol && mod(tstep,10)==0
        plot(x,p,'o'); hold on; plot(xp,Vp*p,'-'); 
        
        plot(xp,Vp*Phi,'--')
        hold off        
        axis([-1 1 -1.5 1.5]); drawnow
    end
    energy(tstep) = sum(sum(wJq.*(Vq*u).^2));
end;

semilogy((1:Nsteps)*dt,energy)
return

%% reverse time

p = 0*p;
u = 0*u;

dt = -dt;
for tstep=1:Nsteps
    for INTRK = 1:5
        id = (tstep-1)*5 + INTRK; % backwards        
        [rhsp rhsu] = WaveRHS1D(p,u,phi,dt,bcOpt,pb(:,id),ub(:,id));
        
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;               
        
        p = p + rk4b(INTRK)*resp;
        u = u + rk4b(INTRK)*resu;        
    end;
       
    if plotSol && (mod(tstep,10)==0 || tstep==Nsteps)
        title('reversal')
        plot(x,p,'o'); hold on; plot(xp,Vp*p,'-'); hold off        
        axis([-1 1 -1.5 1.5]); drawnow
    end
end;

p0 = p; 

% hold on
% plot(xp,Vp*pex(x),'r--')
L2err(1) = sqrt(sum(sum(diag(wq)*(Vq*((p_init-p0).*J)).^2)));   
% title(sprintf('L2 err = %g\n',L2err))

disp(sprintf('L2 err on iter 0 = %0.5g\n',L2err(end)))


%% iterate

% p = initialized to time-reversed condition

p_recon = p0;
p = p0;

for iter = 1:maxiter
        
    p_old = p;
    
    u = zeros(Np,K);
    
    dt = abs(dt);
    for tstep=1:Nsteps
        for INTRK = 1:5
            [rhsp rhsu] = WaveRHS1D(p,u,phi,dt,'abc',[],[]);
            resp = rk4a(INTRK)*resp + dt*rhsp;
            resu = rk4a(INTRK)*resu + dt*rhsu;
            
            p = p + rk4b(INTRK)*resp;
            u = u + rk4b(INTRK)*resu;
            
            id = Nsteps*5 - ((tstep-1)*5 + INTRK) + 1;        
            pb(:,id) = p(vmapB);
            ub(:,id) = u(vmapB);
        end;
        
        if plotSol && mod(tstep,10)==0
            title(sprintf('forward iter %d',iter))
            plot(x,p,'o'); hold on; plot(xp,Vp*p,'-'); 
            axis([-1 1 -1.5 1.5]); drawnow
            hold off
        end
    end;
    
    % go backwards
    if useTraces
        if 1 % compute harmonic extension: does way worse?
            %             R = sparse(N*K+1,(N+1)*K);
            %             for e = 1:K
            %                 ir = (1:N+1) + N*(e-1);
            %                 ic = (1:N+1) + (N+1)*(e-1);
            %                 R(ir,ic) = eye(N+1);
            %             end
            %             Mhat = inv(V*V');
            %             M = kron(spdiag(J(1,:)),Mhat);
            %             Dx = kron(spdiag(rx(1,:)),Dr);
            %             KK = Dx'*M*Dx;
            %             KK = R*KK*R';
            %             pB = zeros(Np,K);
            %             pB(1) = pb(1,1);
            %             pB(end) = pb(2,1);
            %             pCG = diag(1./sum(R,2))*R*pB(:); % nodal averaging for BCs
            %             b = -KK(:,1)*pCG(1) - KK(:,end)*pCG(end);
            %             b(1) = pCG(1);
            %             b(end) = pCG(end);
            %             KK(1,:) = 0; KK(:,1) = 0; KK(1,1) = 1;
            %             KK(end,:) = 0; KK(:,end) = 0; KK(end,end) = 1;
            %
            %             p = reshape(R'*(KK\b),Np,K);
            
            p = (1-x)/2 * pb(1,1) + (1+x)/2 * pb(2,1);
            
            %             keyboard
            %             clf
            %             plot(xp,Vp*p)
            %             return
        else
            p = zeros(size(p)); 
        end
        u = zeros(Np,K);
    end
        
    dt = -dt;
    for tstep=1:Nsteps
        for INTRK = 1:5
            if useTraces
                id = (tstep-1)*5 + INTRK; % backwards
                [rhsp rhsu] = WaveRHS1D(p,u,dt,bcOpt,pb(:,id),ub(:,id));
            else
                [rhsp rhsu] = WaveRHS1D(p,u,phi,dt,'d',[],[]);
            end
            resp = rk4a(INTRK)*resp + dt*rhsp;
            resu = rk4a(INTRK)*resu + dt*rhsu;
            
            p = p + rk4b(INTRK)*resp;
            u = u + rk4b(INTRK)*resu;
        end;
        
        if plotSol && mod(tstep,10)==0            
            plot(x,p,'o'); hold on; plot(xp,Vp*p,'-'); hold off
            axis([-1 1 -1.5 1.5]); drawnow
        end
    end;
    
    if useTraces
        p = p_old - p; % w = (I-A*L)*p0 = p0 - A*L*p0
        p_recon = p_recon + p; % accumulate terms!
    else        
        p = p0 + p;
        p_recon = p;
    end
    
    L2err(iter+1) = sqrt(sum(sum(diag(wq)*(Vq*((p_init-p_recon).*J)).^2)));
    disp(sprintf('L2 err on iter %d = %0.5g\n',iter,L2err(end)))

%     clf
%     plot(x,p,'o');
%     hold on;
%     plot(xp,Vp*p);
%     plot(xp,Vp*p_init,'r--')
    
%     title(sprintf('L2 err = %g\n',L2err(iter)))
end

semilogy(L2err,'o--');hold on

keyboard

function [rhsp rhsu rhsPhi] = WaveRHS1D(p,u,Phi,dt,bcopt,pb,ub)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% fop0 field differences at faces
dp = zeros(Nfp*Nfaces,K);
du = zeros(Nfp*Nfaces,K);
dp(:) = p(vmapP)-p(vmapM);
du(:) = u(vmapP)-u(vmapM);

global tau0
tau = tau0*sign(dt);

if strcmp(bcopt,'abc')
    % ABC
    dp(mapB) = -p(vmapB); 
    du(mapB) = -u(vmapB);
elseif strcmp(bcopt,'d') % read-in boundary values          
    
    dp(mapB) = -2*p(vmapB); 
    du(mapB) = 0;
    
elseif strcmp(bcopt,'rd') % read-in boundary values          
    
    dp(mapB) = 2*(pb-p(vmapB)); 
    du(mapB) = 0;
    
elseif strcmp(bcopt,'rc') % read-in boundary values          
    
    dp(mapB) = (pb-p(vmapB)); 
    du(mapB) = (ub-u(vmapB));
    
end

pflux = (tau*dp - du.*nx);
uflux = (tau*du.*nx - dp).*nx;

% apply characteristic BCs regardless for ABCs
pflux(mapB) = (sign(dt)*dp(mapB) - du(mapB).*nx(mapB));
uflux(mapB) = (sign(dt)*du(mapB).*nx(mapB) - dp(mapB)).*nx(mapB);

% compute right hand sides of the semi-discrete PDE
rhsp = -rx.*(Dr*u) + .5*LIFT*(Fscale.*pflux);
rhsu = -rx.*(Dr*p) + .5*LIFT*(Fscale.*uflux);

global Pq Vq c2 kappa0
kappa = kappa0*sign(dt);
rhsp = Pq*(c2.*(Vq*rhsp) - kappa.*(Vq*p));
rhsp = rhsp - 100*Phi;
rhsPhi = -100*Phi + p; % solve Phi'(t) = -alpha*Phi(t) + p;

return


