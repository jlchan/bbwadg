function L2err = burgers1D_wadg(Nin,K1D)

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;

if nargin==0
    N = 3;
    K1D = 16;
else
    N = Nin;
end
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

vmapP(1) = vmapM(end); % make periodic
vmapP(end) = vmapM(1);

rp = linspace(-1+1e-4,1-1e-4,50)';
Vp = Vandermonde1D(N,rp)/V;
xp=  Vp*x;

global uex Vq Pq

Nq = ceil((3*N-1)/2); % exact integration
% Nq = 2*N+2; % same as overkill integration
Nq = N; % inexact integration
[rq wq] = JacobiGQ(0,0,Nq);
Vq = Vandermonde1D(N,rq)/V;
Pq = V*V'*Vq'*diag(wq);
xq =  Vq*x;

% Set initial conditions
shockSol = 1;
if shockSol
    uex = @(x) 1 + (x < (-1 + .25));
    FinalTime = .75;
    %     FinalTime = 1;
    
else
    % uex = @(x) 1+cos(pi*x);
    uex = @(x) -sin(pi*x);
    FinalTime = .125;
end
u = limit(uex(x));
u = uex(x);


global Wq
%W = @(x) 1./(2+x);
%W = @(x) 1./(1+x.^2);
%W = @(x) exp(-x);
W = @(x) exp(-x.^2);
% W = @(x) 1 + (x > 0);
Wq = W(xq);

% % restore local conservation
[rq2 wq2] = JacobiGQ(0,0,Nq); % exact integration
% [rq2 wq2] = JacobiGQ(0,0,2*N+2); % overkill integration
Vq2 = Vandermonde1D(N,rq2)/V;
Pq2 = V*V'*Vq2'*diag(wq2);
Wq2 = W(Vq2*x);
global Mw
Mw = cell(K,1);
for ee = 1:K
    Mw{ee} = (Pq2*(diag(Wq2(:,ee))*Vq2));
end

global VWL VWR
M = inv(V*V');
invM = V*V';
e = ones(N+1,1);
VWR = zeros(N+1,K);
VWL = zeros(N+1,K);
Ve = zeros(N+1,K);
for ee = 1:K
    Mweight = Vq2'*diag(wq2.*Wq2(:,ee))*Vq2;
    Mwadg = M*((Vq'*diag(wq./Wq(:,ee))*Vq)\M);
    v = (Mwadg-Mweight)*e;
    
    alpha = sign(e'*v)/(e'*v);
    if abs(v'*e) < 1e-16
        alpha = 0;
    end
    v = sqrt(alpha)*v;
    vtilde = Mwadg\v;
    if abs(1+alpha*v'*vtilde) < 1e-14
        keyboard
    end
    
    alphaV = 1/(1 + v'*vtilde);
    VWR(:,ee) = sqrt(alphaV)*M*vtilde;
    VWL(:,ee) = sqrt(alphaV)*vtilde;
    Ve(:,ee) = v;
end

norm(VWL,'fro')
norm(VWR,'fro')


% testing
for ee = 1:K
    Mweight = Vq2'*diag(wq2.*Wq2(:,ee))*Vq2;
    Mcons = M*((Vq'*diag(wq./Wq(:,ee))*Vq)\M) + Ve(:,ee)*Ve(:,ee)'; % conservation fix
    
    invMcons = invM*Vq'*diag(wq./Wq(:,ee))*Vq*invM - VWL(:,ee)*VWR(:,ee)'*invM;
    e1 = norm(invMcons-inv(Mcons),'fro');
    e2 = abs(e'*(Mweight*u(:,ee)) - e'*(Mcons*u(:,ee))); % check conservation holds
    e3 = abs(e'*(Mweight*u(:,ee)) - e'*(invMcons\u(:,ee))); % check inverse works
    if (e1+e2+e3)>1e-12
        keyboard
    end
end
% keyboard
% return

%uexq = fsolve(@(u) u - uex((2+xq).*exp(-u*FinalTime)+2),uex((2+xq).*exp(-uex(xq)*FinalTime)+2));
%uexq = fsolve(@(u) u - uex(tan(atan(xq) - u*FinalTime)),uex(tan(atan(xq) - uex(xq)*FinalTime)));
%uexq = fsolve(@(u) u - uex(-log(exp(-xq)+u*FinalTime)),uex(-log(exp(-xq)+uex(xq)*FinalTime)));
time = FinalTime;
%uexq = fsolve(@(u) u - uex(erfinv((sqrt(pi)*erf(xq)-2*u*time)/sqrt(pi))), uex(erfinv((sqrt(pi)*erf(xq)-2*uex(xq)*time)/sqrt(pi))));
uexq = fsolve(@(u) u - uex(xq-u./W(xq)*time),uex(xq-time));
% plot(xq,uexq)

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K);

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));

dt   = .15*xmin;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;

% outer time step loop
for tstep=1:Nsteps
    
    if shockSol
        % ssp-rk3
        ssprk = [1 .25 2/3];
        utmp = u;
        for i = 1:3
            [rhsu] = BurgersRHS1D(utmp);
            utmp = (1-ssprk(i))*u + ssprk(i)*(utmp + dt*rhsu);
            [utmp Klim alpha] = limit(utmp);
        end
        u = utmp;
    else
        % low storage RK
        for INTRK = 1:5
            rhsu = BurgersRHS1D(u);
            resu = rk4a(INTRK)*resu + dt*rhsu;
            u = u + rk4b(INTRK)*resu;
        end;
    end
    
    time = (tstep-1)*dt;
    
    if nargin==0 && (mod(tstep,5)==0 || tstep==Nsteps)
        
        plot(xp,Vp*u,'b-','linewidth',2)
        hold on;
%         plot(x,u,'o')
        %         plot(xe,VB\u,'o')
        
        if shockSol
            
            xpc = erfinv((sqrt(pi)*erf(xp)-2*1.5*time)/sqrt(pi));
            plot(xp,2 - (xpc > -3/4),'r--','linewidth',2);
            
%                         plot(xp,uex(xp.*W(xp)-1.5*time),'--','linewidth',2)
        end
        hold off
        axis([-1 1 -1 3])
        title(sprintf('Time = %f\n',tstep*dt))
        drawnow
    end
    
end;

if ~shockSol
    err = diag(wq)*((Vq*J).*(Vq*u-uexq).^2);
    L2err = sqrt(sum(err(:)));
    L2err
    if nargin==0
        clf
        plot(xp,Vp*u,'-','linewidth',2)
        hold on;
        plot(x,u,'o')
        
        plot(xq,uexq,'x','linewidth',2);
        hold off
        axis([-1 1 -1 3])
        title(sprintf('Time = %f, error = %g\n',time,L2err))
    end
else
    if nargin==0
        clf
        vv = Vp*u;
        plot(xp(:),vv(:),'b-','linewidth',2,'DisplayName','DG solution')
        hold on;
%         plot(x,u,'.')
        
        xpc = erfinv((sqrt(pi)*erf(xp)-2*1.5*time)/sqrt(pi));
        plot(xp(:),2 - (xpc(:) > -3/4),'r--','linewidth',2,'DisplayName','Exact solution');
%         plot(xp(:),uex(xp(:).*W(xp(:))-1.5*time),'--','linewidth',2,'DisplayName','Exact solution')
        
        hold off
        axis([-1 1 0 2.5])
        title('')
        set(gca,'fontsize',15)
        grid on
        legend show 
    end
    L2err = 0;
end
% keyboard



function [rhsu] = BurgersRHS1D(u)

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
uP(end) = uex(1); % BCs

if 1
    aM = uM;
    du(:) = (uM-uP).*(aM.*nx(:)-(1-alpha)*abs(aM.*nx(:)))/2;
    %     rhsu = -u.*(rx.*(Dr*u)) + LIFT*(Fscale.*(du));
    rhsu = -Pq*((Vq*u).*(Vq*(rx.*(Dr*u)))) + LIFT*(Fscale.*(du));
    
%         global Mw
%         for e=1:K
%             rhsu(:,e) = Mw{e}\rhsu(:,e);
%         end
    
    global Wq
    rhsu = Pq*((Vq*rhsu)./Wq);
    
%         %     locally conservative fix
%         global VWL VWR
%         rhsu = rhsu - VWL.*repmat(sum(VWR.*rhsu,1),N+1,1);
    
else
    %     aM = uM;
    %     du(:) = (uM-uP).*(aM.*nx(:)-(1-alpha)*abs(aM.*nx(:)))/2;
    du(:) = sign(uM).*(uM-uP).*(sign(uM).*nx(:)-abs(sign(uM).*nx(:)))/2;
    %     du(:) = sign(uM).*(uM-uP).*(sign(uM)-abs(sign(uM))).*nx(:)/2;
    rhsu = -(rx.*(Dr*u)) + LIFT*(Fscale.*(du));
    
    global Wq
    rhsu = Pq*((Vq*rhsu).*(Vq*u)./(Wq)); % wadg step
    
end
return



function [u Klim alpha] = limit(u)

[u Klim alpha] = p1limit(u);


function [ulimit ids alpha] = p1limit(u)

Globals1D
global VB

% Compute cell averages
uh = invV*u;
uh(2:Np,:)=0;
uavg = V*uh;
v = uavg(1,:);

% Apply slope limiter as needed.
ulimit = u;
eps0=1.0e-8;
% find end values of each element
ue1 = u(1,:);
ue2 = u(end,:);

% find cell averages
vk = v;
%vkm1 = [v(1),v(1:K-1)]; vkp1 = [v(2:K),v(K)];
vkm1 = [v(K),v(1:K-1)]; vkp1 = [v(2:K),v(1)]; % periodic

% Apply reconstruction to find elements in need of limiting
ve1 = vk - minmod([(vk-ue1);vk-vkm1;vkp1-vk]);
ve2 = vk + minmod([(ue2-vk);vk-vkm1;vkp1-vk]);
ids = find(abs(ve1-ue1)>eps0 | abs(ve2-ue2)>eps0);

if 0
    % convert to bernstein, find K with large TV
    uB = VB\u;
    TV = 0;
    for i = 1:N
        TV = TV + abs(uB(i,:) - uB(i+1,:));
    end
    TV = TV./(N*max(abs(uB(:))));
    ids = find(TV > .5*max(TV));
end

% Check to see if any elements require limiting
if(~isempty(ids))
    % create piecewise linear solution for limiting on specified elements
    uhl = invV*u(:,ids); uhl(3:Np,:)=0; ul = V*uhl;
    % apply slope limiter to selected elements
    ulimit(:,ids) = SlopeLimitLin(ul,x(:,ids),vkm1(ids),vk(ids),vkp1(ids));
end
alpha = 0*ids;

return;

function ulimit = SlopeLimitLin(ul,xl,vm1,v0,vp1)
% function ulimit = SlopeLimitLin(ul,xl,vm1,v0,vp1);
% Purpose: Apply slopelimited on linear function ul(Np,1) on x(Np,1)
%          (vm1,v0,vp1) are cell averages left, center, and right
Globals1D;

% Compute various geometric measures
% ulimit = ul;
h = xl(Np,:)-xl(1,:);
x0 = ones(Np,1)*(xl(1,:) + h/2);
hN = ones(Np,1)*h;

% Limit function
ux = (2./hN).*(Dr*ul);
ulimit = ones(Np,1)*v0+(xl-x0).*(ones(Np,1)*minmod([ux(1,:); (vp1-v0)./h; (v0-vm1)./h]));

return


function mfunc = minmod(v)
% function mfunc = minmodB(v,M,h)
% Purpose: Implement the TVB modified midmod function. v is a vector

Globals1D
M = 50; h = VX(2)-VX(1);
mfunc = v(1,:);
ids = find(abs(mfunc) > M*h.^2);
if(size(ids,2)>0)
    mfunc(ids) = minmod_orig(v(:,ids));
end
return


function mfunc = minmod_orig(v)
% function mfunc = minmod(v)
% Purpose: Implement the midmod function v is a vector
m = size(v,1);
mfunc = zeros(1,size(v,2));
s = sum(sign(v),1)/m;
ids = find(abs(s)==1);
if(~isempty(ids))
    mfunc(ids) = s(ids).*min(abs(v(:,ids)),[],1);
end
return;