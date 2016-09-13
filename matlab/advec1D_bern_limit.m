function advec1D_bern_limit

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

global VB DB VBe W
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
u = limit(uex(x));
% u = uex(x);

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
        [rhsu] = AdvecRHS1D(utmp);
        plot(xp,Vp*u,'-')
        utmp = (1-ssprk(i))*u + ssprk(i)*(utmp + dt*rhsu);
        [utmp Klim alpha] = limit(utmp);
    end    
    u = utmp;

    % Increment time
    time = time+dt;
    
    if mod(tstep,5)==0
        
        plot(xp,Vp*u,'-')
        %         plot(x,VB*uex(xe-time),'.-')
        hold on;
        plot(x,u,'o')
%         plot(xe,VB\u,'o')
        
        plot(xp,uex(xp-time),'--');
        if (any(alpha>1e-8))
            plot(x(:,Klim),alpha,'x')
        end
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

function [rhsu] = AdvecRHS1D(u)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form field differences at faces
alpha = 0;
du = zeros(Nfp*Nfaces,K);
du(:) = (u(vmapM)-u(vmapP)).*(nx(:)-(1-alpha)*abs(nx(:)))/2;

% compute right hand sides of the semi-discrete PDE
rhsu = -rx.*(Dr*u) + LIFT*(Fscale.*(du));
return


function [u Klim alpha] = limit(u)

[u Klim alpha] = BBlimit(u);
% [u Klim alpha] = p1limit(u);

function [u Klim alpha] = BBlimit(u)

Globals1D
global VB DB VBe W

% convert to bernstein
uB = VB\u;

% find K with large TV
TV = 0;
for i = 1:N
    TV = TV + abs(uB(i,:) - uB(i+1,:));
end
TV = TV./(N*max(abs(uB(:))));
Klim = find(TV > .75*max(TV));
uB = uB(:,Klim);

% % local conservation
% uavg = repmat(sum(uB,1)/(N+1),N+1,1);

if 1
    [u1] = p1limit(u); u1 = u1(:,Klim);
else
    % limit P1 part
    u1 = V\u(:,Klim); uavg = u1(1,:); u1(3:end,:) = 0; u1 = V*u1;
    umax = 2; umin = 1; pmax = max(u1,[],1); pmin = min(u1,[],1);
    theta = min(abs([(umax - uavg)./(pmax-uavg);(umin-uavg)./(pmin-uavg)]),[],1);
    theta = min([theta;ones(size(theta))]);
    uavg = repmat(uavg,N+1,1);
    
    % limited P1 sol
    u1 = (u1 - uavg)*diag(theta) + uavg;
end

% u1 = 0*u1;

% remove P1 part
uB = uB - VB\u1;

if 1
%     uBe = VBe*uB; % BB approx
    uBe = W*uB; % locally conservative BB
    err = max(abs(uB-uBe));
    
    alpha = zeros(size(err));
    if err > max(abs(uB))/sqrt(N) % bound variation?
        alpha = err/sqrt(N);
    end
%     alpha(:) = 1;
    alpha = min(1,alpha);
    uB = uB*diag(1-alpha) + uBe*diag(alpha);
else
    uB = 0*uB;
    alpha = zeros(K,1);
end
u(:,Klim) = u1 + VB*uB; % convert back


function [ulimit ids alpha] = p1limit(u)

Globals1D

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
ve1 = vk - minmodB([(vk-ue1);vk-vkm1;vkp1-vk]);
ve2 = vk + minmodB([(ue2-vk);vk-vkm1;vkp1-vk]);
ids = find(abs(ve1-ue1)>eps0 | abs(ve2-ue2)>eps0);

if 0
    global VB
    
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
alpha = ones(size(ids))*.5;

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
ulimit = ones(Np,1)*v0+(xl-x0).*(ones(Np,1)*minmodB([ux(1,:); (vp1-v0)./(h/2); (v0-vm1)./(h/2)]));

return


function mfunc = minmodB(v)
% function mfunc = minmodB(v,M,h)
% Purpose: Implement the TVB modified midmod function. v is a vector

Globals1D
M = 50; h = VX(2)-VX(1);
mfunc = v(1,:);
ids = find(abs(mfunc) > M*h.^2);
if(size(ids,2)>0)
    mfunc(ids) = minmod(v(:,ids));
end
return


function mfunc = minmod(v)
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