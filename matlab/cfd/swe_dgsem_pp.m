clear

global N K

N = 5;
K = 32; % number of elements
CFL = .125; % should be O(1) - too large, soln blows up
FinalTime = 4;
tau = 1;

%% nodal DG

global Vf PN DNrskew L mapP nx J

[r w] = JacobiGL(0,0,N);
rq = r; wq = w;

V = Vandermonde1D(N,r);
Dr = GradVandermonde1D(N,r)/V;
Vq = Vandermonde1D(N,rq)/V;
Vf = Vandermonde1D(N,[-1;1])/V;

M = Vq'*diag(wq)*Vq;
L = M\(Vf');
Qr = M*Dr;
Qrskew = (Qr-Qr');

QNr = [Qrskew Vf'*diag([-1;1])
    -diag([-1;1])*Vf zeros(2)];
WN = diag([wq;1;1]);
PN = M\[Vq' Vf']*WN;
DNrskew = WN\QNr;

DNrskew = M\(Qrskew); % for LGL only! 

VX = linspace(-1,1,K+1)';
h = VX(2)-VX(1);
x = zeros(N+1,K);
for e = 1:K
    x(:,e) = h*(1+r)/2 + VX(e);
end
J = h/2;

nx = repmat([-1;1],1,K);

mapP = zeros(2,K);
for e = 2:K-1
    mapP(:,e) = [-2;1] + 2*e;
end
mapP(:,1) = [2*K; 3];
mapP(:,end) = [2*K-2; 1];

%% make fv sub-cells

global dx mapL mapR mapPfv nxfv

rf = [-1;cumsum(w)-1];
dx = h*repmat(w,1,K)/sum(w);

% map for L/R cell centers
mapL = reshape([(N+1)*K 1:(N+1)*K-1],N+1,K);
mapR = reshape([2:(N+1)*K 1],N+1,K);

mapPfv = zeros(2,(N+1)*K);
for e = 2:(N+1)*K-1
    mapPfv(:,e) = [-2;1] + 2*e;
end
mapPfv(:,1) = [2*(N+1)*K; 3];
mapPfv(:,end) = [2*(N+1)*K-2; 1];

nxfv = repmat([-1;1],1,(N+1)*K);

%%

global V1 V2 U1 U2 fSn1 fSn2 fS1 fS2
global fSnfv1 fSnfv2
global fSn fS

% advection
fS = @(uL,uR) .5*(uL + uR);
fD = @(uL,uR) .5*(uR-uL);

% % Burgers
% fS = @(uL,uR) (uL.^2 + uL.*uR + uR.^2)/6;
% fD = @(uL,uR) .5*abs(.5*(uL+uR)).*(uR-uL);

% shallow water
g = 1; % assume b = const for now
Sfun = @(h,u) .5*(h.*u.^2 + g*h.^2);
V1 = @(h,u) g*h - .5*u.^2;
V2 = @(h,u) u;
U1 = @(v1,v2) (v1+.5*v2.^2)/g;
U2 = @(v1,v2) v2;
fS1 = @(hL,hR,uL,uR) (hL.*uL + hR.*uR)/2;
fS2 = @(hL,hR,uL,uR) .5*(hL.*uL + hR.*uR).*.5.*(uL+uR) + .5*g.*(hL.*hR);

% Lfc = @(hL,hR,uL,uR) abs(.5*(uL+uR))+sqrt(g*.5*(hL+hR)); % |u| + sqrt(gh)
Lf = @(h,u) abs(u)+sqrt(g*h);
Lfc = @(hL,hR,uL,uR) max(Lf(hL,uL),Lf(hR,uR)); % stronger LF
fD1 = @(hL,hR,uL,uR) .5*Lfc(hL,hR,uL,uR).*(hR-hL);
fD2 = @(hL,hR,uL,uR) .5*Lfc(hL,hR,uL,uR).*(hR.*uR-hL.*uL);

fSn1 = @(hL,hR,uL,uR,nx) fS1(hL,hR,uL,uR).*nx - tau*fD1(hL,hR,uL,uR);
fSn2 = @(hL,hR,uL,uR,nx) fS2(hL,hR,uL,uR).*nx - tau*fD2(hL,hR,uL,uR);

% interface flux
fSn = @(uL,uR,nx) fS(uL,uR).*nx - tau*fD(uL,uR); %.5*(uL.^b+uR.^b).*nx - tau*.5*(uR-uL).*abs(.5*(uR+uL)).^(b-1);

%% initial solution

h0 = @(x) rand(N+1,K)+1e-8; %
u0 = @(x) randn(N+1,K);

h0 = @(x) 1e-1 + exp(-50*x.^2); % gaussian initial condition
u0 = @(x) 0*exp(-50*x.^2); % gaussian initial condition

h0 = @(x) .5*(x < 0) + .5 + 0*x;

h = h0(x);
hu = h.*u0(x);

% h = [    0.0010    0.0010    0.0056    0.0007    0.0007    0.0050    0.5414    1.0188    1.0056    1.0369    0.4493    0.0044    0.0031    0.0050    0.0011    0.0014
%     0.0011    0.0010    0.0003    0.0013    0.0013    0.0787    0.7347    0.9871    1.0019    0.9585    0.2944    0.0015    0.0006    0.0005    0.0010    0.0009
%     0.0009    0.0010    0.0005    0.0006    0.0015    0.2944    0.9585    1.0019    0.9871    0.7347    0.0787    0.0013    0.0013    0.0003    0.0010    0.0011
%     0.0014    0.0011    0.0050    0.0031    0.0044    0.4493    1.0369    1.0056    1.0188    0.5414    0.0050    0.0007    0.0007    0.0056    0.0010    0.0010];
% hu = [    0.0000   -0.0000    0.0024   -0.0001   -0.0345   -0.0033   -0.2869    0.0200   -0.0135   -0.0390    0.3064   -0.0034    0.0002    0.0019   -0.0003    0.0000
%    -0.0000    0.0000   -0.0007    0.0000    0.0154   -0.0447   -0.2097   -0.0167   -0.0039    0.0420    0.2253    0.0075    0.0003   -0.0006    0.0001   -0.0000
%     0.0000   -0.0001    0.0006   -0.0003   -0.0075   -0.2253   -0.0420    0.0039    0.0167    0.2097    0.0447   -0.0154   -0.0000    0.0007   -0.0000    0.0000
%    -0.0000    0.0003   -0.0019   -0.0002    0.0034   -0.3064    0.0390    0.0135   -0.0200    0.2869    0.0033    0.0345    0.0001   -0.0024    0.0000   -0.0000];
%
% S = h.^2 + h.*(hu./h).^2;

%% solver

rk4a = [            0.0 ...
    -567301805773.0/1357537059087.0 ...
    -2404267990393.0/2016746695238.0 ...
    -3550918686646.0/2091501179385.0  ...
    -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
    5161836677717.0/13612068292357.0 ...
    1720146321549.0/2090206949498.0  ...
    3134564353537.0/4481467310338.0  ...
    2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
    1432997174477.0/9575080441755.0 ...
    2526269341429.0/6820363962896.0 ...
    2006345519317.0/3224310063776.0 ...
    2802321613138.0/2924317926251.0];

ssprkb = [1 .25 2/3];

% CN = (N+1)^2;
lam = abs(hu./h) + sqrt(g*h); %abs(u(:));
dt = CFL*min(dx(:))/max(lam(:)); % set time-step
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

figure(1)
resfv1 = zeros(N+1,K);
resfv2 = zeros(N+1,K);
res1 = zeros(N+1,K);
res2 = zeros(N+1,K);

for i = 1:Nsteps
    
    hh = h;
    hhuu = hu;
    for INTRK = 1:3
        
        b = ssprkb(INTRK);

        [rhs1, rhs2] = rhs_DG(hh,hhuu);
        hDG = (1-b)*h + b*(hh + dt*rhs1);
        huDG = (1-b)*hu + b*(hhuu + dt*rhs2);
        
        % cell averages instead of FV
        hfv = repmat(.5*(w'*hDG),N+1,1);
        hufv = repmat(.5*(w'*huDG),N+1,1);        
        hmin = repmat(min(hDG),N+1,1);
        a = max(min(1,hfv./(1e-10+hfv-hmin)),0);       

        if min(hfv(:))<0
            plot(x,h,'bo-')
            hold on
            plot(x,hfv,'r--') 
        end

        
        hh = hfv + a.*(hDG-hfv);
        hhuu = hufv + a.*(huDG-hufv);
        

        %         if i==21 && INTRK==2
        %             clf
%                     plot(x,h,'bo--')
%                     hold on
%                     plot(x,hfv,'r^')
        %             title(sprintf('timestep %d\n',i))
        %             keyboard
        %         end
    end
    h = hh;
    hu = hhuu;
    
    
    
    %     for INTRK = 1:5
    %
    %         [rhsfv1 rhsfv2] = rhs_fv(h,hu);
    %         resfv1 = rk4a(INTRK)*resfv1 + dt*rhsfv1;
    %         resfv2 = rk4a(INTRK)*resfv2 + dt*rhsfv2;
    %         hfv = h + rk4b(INTRK)*resfv1;
    %         hufv = hu + rk4b(INTRK)*resfv2;
    %
    %         [rhs1, rhs2] = rhs_DG(h,hu);
    %         res1 = rk4a(INTRK)*res1 + dt*rhs1;
    %         res2 = rk4a(INTRK)*res2 + dt*rhs2;
    %         h = h + rk4b(INTRK)*res1;
    %         hu = hu + rk4b(INTRK)*res2;
    %
    %         hDG = ssprku(INTRK)*h + (1-ssprku(INTRK))*hh + ssprkb(INTRK)*dt*rhs1;
    %         hhuu = ssprku(INTRK)*hu + (1-ssprku(INTRK))*hhuu + ssprkb(INTRK)*dt*rhs2;
    %
    %         a = max(min(1,min(hfv)./(min(hfv)-min(h))),0);
    %         a(min(h) > 1e-12) = 1;
    % %         if norm(a-1)>1e-8
    % %             keyboard
    % %         end
    %         a = repmat(a,N+1,1);
    % %         a = 0;
    %         h = hfv + a.*(h-hfv);
    %         hu = hufv + a.*(hu-hufv);
    %     end
    
    if i==Nsteps || (mod(i,25)==0)
        %         plot(x,u,'bo--')
        plot(x,h,'bo--')
        hold on
        plot(x,hu,'r^')
        plot(x,a,'k*')
        hold off
        title(sprintf('min h = %g\n',min(h(:))))
        axis([-1 1 -.1 1.5])
        drawnow
    end
end

% figure(2)
% hold on
% plot(x,u,'x--')

function [rhs1 rhs2] = rhs_DG(h,hu)

global N K
global Vf PN DNrskew L mapP nx J
global V1 V2 U1 U2 fSn1 fSn2 fS1 fS2

u = hu./h;
hM = h([1;N+1],:); hP = hM(mapP);
uM = hu([1;N+1],:)./hM; uP = uM(mapP);
flux1 = fSn1(hM,hP,uM,uP,nx);
flux2 = fSn2(hM,hP,uM,uP,nx);

rhs1 = zeros(N+1,K);
rhs2 = zeros(N+1,K);
for e = 1:K
    [hx hy] = meshgrid(h(:,e));
    [ux uy] = meshgrid(u(:,e));
    rhs1(:,e) = -(sum(DNrskew.*fS1(hx,hy,ux,uy),2) + L*flux1(:,e))./J;
    rhs2(:,e) = -(sum(DNrskew.*fS2(hx,hy,ux,uy),2) + L*flux2(:,e))./J;
end


end
