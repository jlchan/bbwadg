function L2err = burgers1D_wadg(Nin,K1D)

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;

if nargin==0
    N = 4;
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

global VB DB VBe W
re = linspace(-1,1,N+1)';
[VB VBr] = bern_basis_1D(N,r);
DB = VB\VBr;
DB(abs(DB)<1e-8) = 0;
VBe = bern_basis_1D(N,re);
xe = Vandermonde1D(N,re)/V * x;

W = get_BB_smoother(N);
% W = bern_basis_1D(N,re);

global uex Vq Pq

Nq = ceil((3*N-1)/2); % exact integration
% Nq = 2*N+2; % same as overkill integration
Nq = N; % inexact integration
[rq wq] = JacobiGQ(0,0,Nq);
Vq = Vandermonde1D(N,rq)/V;
Pq = V*V'*Vq'*diag(wq);
xq =  Vq*x;

% Set initial conditions
uex = @(x) 1 + (x < (-1 + 1/K1D));
% uex = @(x) -sin(pi*x);
% uex = @(x) -1 - (x < 0);
% uex = @(x) 1-cos(pi*x);
% uex = @(x) exp(-100*x.^2);
FinalTime = 1.2;
u = limit(uex(x));
% u = uex(x);

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K);

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));

dt   = .15*xmin;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;

% outer time step loop
for tstep=1:Nsteps
    
    % ssp-rk3
    ssprk = [1 .25 2/3];
    utmp = u;
    for i = 1:3
        [rhsu] = BurgersRHS1D(utmp);
        utmp = (1-ssprk(i))*u + ssprk(i)*(utmp + dt*rhsu);
        [utmp Klim alpha] = limit(utmp);
    end
    u = utmp;
        
    time = (tstep-1)*dt;
    
    if nargin==0 && (mod(tstep,5)==0 || tstep==Nsteps)
        
        plot(xp,Vp*u,'b-','linewidth',2)
        hold on;
        %         plot(x,u,'o')
        %         plot(xe,VB\u,'o')
        plot(VX,ones(size(VX)),'x')
%         plot(xp,uex(xp-1.5*time),'--','linewidth',2)
        hold off
        axis([-1 1 -1 3])
        title(sprintf('Time = %f\n',tstep*dt))
        drawnow
    end
    
end;



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
% uP(end) = uex(1); % BCs

opt = 1;
if opt==1
    aM = uM;
    du(:) = (uM-uP).*(aM.*nx(:)-(1-alpha)*abs(aM.*nx(:)))/2;
        rhsu = -u.*(rx.*(Dr*u)) + LIFT*(Fscale.*(du));
%     rhsu = -Pq*((Vq*u).*(Vq*(rx.*(Dr*u)))) + LIFT*(Fscale.*(du));    
elseif opt==2 % advec only
    a = 1;
    aM = a;
    du(:) = (uM-uP).*(aM.*nx(:)-(1-alpha)*abs(aM.*nx(:)))/2;
    rhsu = -a.*(rx.*(Dr*u)) + LIFT*(Fscale.*(du));
else
    du(:) = sign(uM).*(uM-uP).*(sign(uM).*nx(:)-abs(sign(uM).*nx(:)))/2;
    rhsu = -(rx.*(Dr*u)) + LIFT*(Fscale.*(du));
    
    rhsu = Pq*((Vq*rhsu).*(Vq*u)); % wadg step
end

return


function [u Klim alpha] = limit(u)

% [u Klim alpha] = p1limit(u);
[u Klim alpha] = BBlimit(u);

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
Klim = find(TV > .5*max(TV));
% Klim = 1:K;
uB = uB(:,Klim);

opt=1;
if opt==1 % default to minmod for p1
    [u1] = p1limit(u); u1 = u1(:,Klim); 
elseif opt==2
    % limit P1 part
    u1 = V\u(:,Klim); uavg = u1(1,:); u1(3:end,:) = 0; u1 = V*u1;
    umax = 2; umin = 1; pmax = max(u1,[],1); pmin = min(u1,[],1);
    theta = min(abs([(umax - uavg)./(pmax-uavg);(umin-uavg)./(pmin-uavg)]),[],1);
    theta = min([theta;ones(size(theta))]);
    uavg = repmat(uavg,N+1,1);
    
    % limited P1 sol
    u1 = (u1 - uavg)*diag(theta) + uavg;
else
    u1 = 0*uB; % ignore P1 part?
end

% remove P1 part
uB = uB - VB\u1;

if 1
    % uBe = VBe*uB; % BB approx
    uBe = W*uB;
    % err = max(abs(uB - uBe))*diag(1./max(abs(DB*DB*uB)));
    err = max(abs(uB-uBe));
    
    alpha = zeros(size(err));
    if err > max(abs(uB))/sqrt(N) % bound variation?
        alpha = err/sqrt(N);
    end
    alpha = min(1,alpha);    
    uB = uB*diag(1-alpha) + uBe*diag(alpha);
else
    uB = 0*uB;
    alpha = zeros(K,1);
end
u(:,Klim) = u1 + VB*uB; % convert back to nodal



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