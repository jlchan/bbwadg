% wave

function Wave1D_SEM

% Driver script for solving the 1D advection equations
Globals1D;

% Order of polymomials used for approximation
N = 8;

CG = 0;
% Solve Problem
FinalTime = 8;

% Generate simple mesh
K1D = 8;
% K1D = 64;
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

% s = [ones(N-1,1); .5; 0];
% a = .5; kc = 1;
% s = ones(N+1,1); s(kc:N+1) = 1-a*(((kc:N+1) - kc)./(N-kc)).^2;
% plot(0:N,s);return
% F = V*diag(s)/V;

re = linspace(-1,1,N+1)';
Ve = Vandermonde1D(N,re)/V;
xe = Ve*x;
VB = bern_basis_1D(N,r);
VBe = bern_basis_1D(N,re);
F = VB*Ve;
% F = get_BB_P1smoother(N);

TV = TV1D(VB\p) + TV1D(VB\u);
[~,id] = max(TV);
ids = id-2:id+1; % refine in neighborhood

plim = p;
plim(:,ids) = F*plim(:,ids);
% [plim ids] = p1limit(plim);

figure
% plot(x(:),p(:),'o','linewidth',2,'markersize',8,'MarkerFaceColor','g');
hold on
pp = Vp*p;
plot(xp(:),pp(:),'-.','linewidth',2);
plot(xp(:),pex(xp(:)),'-','linewidth',2);
axis([-1 1 -.5 1.5])
hold off
grid on
legend('DG','Exact solution')
set(gca,'fontsize',15)

print(gcf,'-dpng','DGosc.png')

figure
ppl = Vp*plim;
hold on
plot(xp(:),ppl(:),'-.','linewidth',2);
plot(xp(:),pex(xp(:)),'-','linewidth',2);
axis([-1 1 -.5 1.5])
hold off
grid on
legend('Filtered DG','Exact solution')
print(gcf,'-dpng','DGfilter.png')
%legend('FEM approximation','Exact solution')
set(gca,'fontsize',15)
% print(gcf,'-dpng','CG.png')
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


function [ulimit ids] = p1limit(u)

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

% Check to see if any elements require limiting
if(~isempty(ids))
    % create piecewise linear solution for limiting on specified elements
    uhl = invV*u(:,ids); uhl(3:Np,:)=0; ul = V*uhl;
    % apply slope limiter to selected elements
    ulimit(:,ids) = SlopeLimitLin(ul,x(:,ids),vkm1(ids),vk(ids),vkp1(ids));
end

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
M = 40; h = VX(2)-VX(1);
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

