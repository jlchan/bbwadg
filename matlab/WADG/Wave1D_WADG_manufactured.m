clear

Globals1D;
global Pq Vq wq M

% Order of polymomials used for approximation
N = 7;
MM = 7;

FinalTime = .7;
CFL = .125; 

% Kvec = 16;
Kvec = [2:10];
% Kvec = [4 8 16];% 32]; % 64];% 128];

sk = 1;
for K1D = Kvec
    
    K1D
    [Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
    
    % Initialize solver and construct grid and metric
    StartUp1D;    
        
    vmapP(1) = vmapM(end); % hack for periodic
    vmapP(end) = vmapM(1); % hack for periodic
    
    [rq wq] = JacobiGQ(0,0,N+1);
    Vq = Vandermonde1D(N,rq)/V;
    M = Vq'*diag(wq)*Vq;
    Pq = M\(Vq'*diag(wq));
    xq = Vq*x;    
    
    [rq2 wq2] = JacobiGQ(0,0,N+2);
    Vq2 = Vandermonde1D(N,rq2)/V;
    xq2 = Vq2*x;    
    wJq = diag(wq2)*(Vq2*J);
    
    % Set initial conditions
%     p = exp(-50*x.^2);
    p = cos(pi*x);
    u = zeros(Np,K);
    
    % BBwadg version
    p2 = p;
    u2 = u;
   
    c2 = @(x) 1 + .25*sin(.3+pi*x);
%     c2 = @(x) exp(sin(1+2*pi*x)); 
    cq = c2(xq);    
    cq2 = Vandermonde1D(MM,rq)*Vandermonde1D(MM,rq)'*(wq.*cq);
    
    global f
    f = @(t) zeros(Np,K);
    
%     rM = JacobiGL(0,0,MM);
%     VNM = Vandermonde1D(N,rM)/V;
%     VMN = Vandermonde1D(MM,r)/Vandermonde1D(MM,rM);
%     cq2 = Vq*VMN*c2(VNM*x);
    
    %% Solve Problem
    
    time = 0;
    
    ssprkb = [1 .25 2/3];
    
    % Runge-Kutta residual storage
    resp = zeros(Np,K);
    resu = zeros(Np,K);
    resp2 = zeros(Np,K);
    resu2 = zeros(Np,K);    
    
    % compute time step size
    xmin = min(abs(x(1,:)-x(2,:)));
    dt  = CFL*xmin;
    Nsteps = ceil(FinalTime/dt);
    dt = FinalTime/Nsteps;
    
    % outer time step loop
    for tstep=1:Nsteps
        
        if 0 % euler or ssprk3
             
            [rhsp rhsu] = WaveRHS1D(p,u,cq);
            [rhsp2 rhsu2] = WaveRHS1D(p2,u2,cq2);
                        
            p = p + dt*rhsp;
            u = u + dt*rhsu;            
            p2 = p2 + dt*rhsp2;
            u2 = u2 + dt*rhsu2;                        
            
%             pp = p; uu = u;
%             pp2 = p2; uu2 = u2;
%             for INTRK = 1:3                
%                 [rhsp rhsu] = WaveRHS1D(pp,uu,cq);
%                 [rhsp2 rhsu2] = WaveRHS1D(pp2,uu2,cq2);
%                                
%                 b = ssprkb(INTRK);
%                 pp = (1-b)*p + b*(pp + dt*rhsp);
%                 uu = (1-b)*u + b*(uu + dt*rhsu);
%                 pp2 = (1-b)*p2 + b*(pp2 + dt*rhsp2);
%                 uu2 = (1-b)*u2 + b*(uu2 + dt*rhsu2);
%             end
%             p = pp; u = uu;
%             p2 = pp2; u2 = uu2;
        
        else % rk4
            
            for INTRK=1:5
                
                timeloc = time + rk4c(INTRK)*dt;
                
                [rhsp rhsu] = WaveRHS1D(p,u,cq,timeloc);
                [rhsp2 rhsu2] = WaveRHS1D(p2,u2,cq2,timeloc);  
                                                
                resp = rk4a(INTRK)*resp + dt*rhsp;
                resu = rk4a(INTRK)*resu + dt*rhsu;
                p = p + rk4b(INTRK)*resp;
                u = u + rk4b(INTRK)*resu;
                
                resp2 = rk4a(INTRK)*resp2 + dt*rhsp2;
                resu2 = rk4a(INTRK)*resu2 + dt*rhsu2;                                
                p2 = p2 + rk4b(INTRK)*resp2;
                u2 = u2 + rk4b(INTRK)*resu2;
            end;                        
            
        end
        
        
        % Increment time
        time = time+dt;
        
%         pdiff(tstep) = sqrt(sum(sum(wJq.*(Vq2*p-Vq2*p2).^2)));        
        
        if length(Kvec)==1 && mod(tstep,10)==0
            clf
            plot(x,p,'--');
            hold on
            plot(x,p2,'o--');
            axis([-1,1,-2,2])
            drawnow
        end
                
    end
    
    if length(Kvec)==1
        return
    end
    
%     figure(2)
%     plot(xq2,Vq2*(p-p2),'-')
%     return;
    
    L2diff(sk) = sqrt(sum(sum(wJq.*(Vq2*p-Vq2*p2).^2)));
%     L2diff(sk) = (sum(sum(wJq.*(cq./cq2-1).^2)));
%     L2diff(sk) = abs(sum(sum(wJq.*(cq-cq2).*(Vq*rhsp2).*(Vq*(p-p2)))));
    
    sk = sk + 1;
end

%%
figure(2)
h = 2./Kvec;
loglog(h,L2diff,'o--')
hold on
r = min(MM+3,N+1);
% loglog(h,5e-2*h.^r,'k--')
% loglog(h,1e-2*h.^(r-1),'k-.')
rates = diff(log(L2diff))./diff(log(h));
% rates = rates(end-1:end);

title(sprintf('est rate = %f\n',mean(rates)))

scale = L2diff(1)/h(1).^mean(rates);
loglog(h,scale*h.^mean(rates),'k--')

%%
function [rhsp rhsu] = WaveRHS1D(p,u,cq,t)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;
global Pq Vq wq M

% form field differences at faces
dp = zeros(Nfp*Nfaces,K);
du = zeros(Nfp*Nfaces,K);
dp(:) = p(vmapP)-p(vmapM);
du(:) = u(vmapP)-u(vmapM);

if 1
    dp(mapB) = 0;
    du(mapB) = -2*u(vmapB);
end

tau = 0;
pflux = .5*(tau*dp - du.*nx);
uflux = .5*(tau*du.*nx - dp).*nx;

% compute right hand sides of the semi-discrete PDE
global f
rhsp = -rx.*(Dr*u) + LIFT*(Fscale.*pflux) + f(t);
rhsu = -rx.*(Dr*p) + LIFT*(Fscale.*uflux);

rhsp = Pq*(cq.*(Vq*rhsp));
% for e = 1:K 
%     Mw = Vq'*diag(wq./cq(:,e))*Vq;
%     Pw = Mw\M;
%     rhsp(:,e) = Pw*rhsp(:,e);
% end

end

