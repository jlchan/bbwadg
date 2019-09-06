clear

global N K

N = 3;
K = 32; % number of elements
CFL = .25; % should be O(1) - too large, soln blows up
FinalTime = 1;
tau = 1;

%% nodal DG

global Vf PN DNrskew Drw L mapP nx J

[r w] = JacobiGL(0,0,N);
rq = r; wq = w;
% [rq wq] = JacobiGQ(0,0,N+1);

V = Vandermonde1D(N,r);
Dr = GradVandermonde1D(N,r)/V;
Vq = Vandermonde1D(N,rq)/V;
Vf = Vandermonde1D(N,[-1;1])/V;

M = Vq'*diag(wq)*Vq;
% Pq = M\(Vq'*diag(wq));
L = M\(Vf');
Qr = M*Dr;
Qrskew = (Qr-Qr');
QNr = [Qrskew Vf'*diag([-1;1])
    -diag([-1;1])*Vf zeros(2)];
B = [zeros(N+1) 0*Vf';
    0*Vf diag([-1;1])];
WN = diag([wq;1;1]);
PN = M\[Vq' Vf']*WN;

Drw = M\(Dr'*M);

DNrskew = WN\QNr;

% rN = [rq;-1;1];
% f = @(x) exp(x);
% plot(rq,M\[Vq' Vf']*(.5*DNrskew+.5*B)*f(rN),'o--')
% hold on
% plot(rq,f(rq),'x')

%%
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
fS2 = @(hL,hR,uL,uR) .5*(hL.*uL + hR.*uR).*.5.*(uL+uR) + .5*g.(hL.*hR);

% Lfc = @(hL,hR,uL,uR) abs(.5*(uL+uR))+sqrt(g*.5*(hL+hR)); % |u| + sqrt(gh)
Lf = @(h,u) abs(u)+sqrt(g*h);
Lfc = @(hL,hR,uL,uR) max(Lf(hL,uL),Lf(hR,uR)); % stronger LF
fD1 = @(hL,hR,uL,uR) .5*Lfc(hL,hR,uL,uR).*(hR-hL);
fD2 = @(hL,hR,uL,uR) .5*Lfc(hL,hR,uL,uR).*(hR.*uR-hL.*uL);

fSn1 = @(hL,hR,uL,uR,nx) fS1(hL,hR,uL,uR).*nx - tau*fD1(hL,hR,uL,uR);
fSn2 = @(hL,hR,uL,uR,nx) fS2(hL,hR,uL,uR).*nx - tau*fD2(hL,hR,uL,uR);

%% initial solution

h0 = @(x) rand(N+1,K)+1e-8; %
u0 = @(x) randn(N+1,K);

h0 = @(x) .01 + exp(-50*x.^2); % gaussian initial condition
% u0 = @(x) exp(-50*x.^2); % gaussian initial condition

% a = .1;
% h0 = @(x) (1-a)*(x < 0) + a*ones(size(x));
% h0 = @(x) 1*(x < 0) + ones(size(x));
u0 = @(x) zeros(size(x)); % gaussian initial condition

% h0 = @(x) ones(size(x));
% u0 = @(x) (x > 0)-.5;


h = h0(x);
hu = h.*u0(x);
% plot(x,hu,'o')
% return
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

% CN = (N+1)^2;
lam = abs(hu./h) + sqrt(g*h); %abs(u(:));
dt = CFL*min(dx(:))/max(lam(:)); % set time-step
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

wJq = diag(wq)*repmat(J,N+1,K);

figure(1)
resfv1 = zeros(N+1,K);
resfv2 = zeros(N+1,K);
res1 = zeros(N+1,K);
res2 = zeros(N+1,K);

global a ep1 ep2
for i = 1:Nsteps      
    
    for INTRK = 1:5
        
        a = ones(N+1,K);
        c = V\h;
        sN = log(max(sum(c(N+1,:).^2,1)./sum(c.^2,1),sum(c(N,:).^2,1)./sum(c(1:N,:).^2,1)));
        ep1 = 2*J/N; % h/p
        ep2 = 2*J/N; % h/p
        s0 = log(1/N.^4);
        kappa = 10;
        sfun = .5*(1+sin(pi*(sN-s0)/(2*kappa)));
        a(:,sN < s0-kappa) = 0;
        a(:,(sN > s0-kappa) & (sN < s0 + kappa)) = repmat(sfun((sN > s0-kappa) & (sN < s0 + kappa)),N+1,1);
        a(:,(sN > s0+kappa)) = 1;
        a(:,min(h)>.1) = 0;
%         a(:) = 0;
        
        [rhs1, rhs2] = rhs_DG(h,hu);
        
%         sum(sum(wJq.*(V1(h,hu./h).*rhs1 + V2(h,hu./h).*rhs2)))
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        h = h + rk4b(INTRK)*res1;
        hu = hu + rk4b(INTRK)*res2;
        
    end
    
    if i==Nsteps || (mod(i,5)==0)
        %         plot(x,u,'bo--')
        plot(x,h,'bo-')
        hold on
        plot(x,hu,'r^-')
        plot(x,a,'ko-')
        %         plot(x,Lf,'ks')
        hold off
        title(sprintf('min h = %g at time %f\n',min(h(:)),dt*i))
        %         axis([-1 1 -2 2.5])
        drawnow
    end
end

% figure(2)
% hold on
% plot(x,u,'x--')

function [rhs1 rhs2] = rhs_DG(h,hu)

global N K
global Vf Pq PN DNrskew Drw L mapP rxJ nx J
global V1 V2 U1 U2 fSn1 fSn2 fS1 fS2

u = hu./h;
v1M = Vf*V1(h,u);
v2M = Vf*V2(h,u);
hM = U1(v1M,v2M); hP = hM(mapP);
uM = U2(v1M,v2M); uP = uM(mapP);
flux1 = fSn1(hM,hP,uM,uP,nx);
flux2 = fSn2(hM,hP,uM,uP,nx);

hN = [h;hM];
uN = [u;uM];
rhs1 = zeros(N+1,K);
rhs2 = zeros(N+1,K);
for e = 1:K
    [hx hy] = meshgrid(hN(:,e));
    [ux uy] = meshgrid(uN(:,e));
    rhs1(:,e) = -(PN*sum(DNrskew.*fS1(hx,hy,ux,uy),2) + L*flux1(:,e))./J;
    rhs2(:,e) = -(PN*sum(DNrskew.*fS2(hx,hy,ux,uy),2) + L*flux2(:,e))./J;
end


% viscous stresses:
hM = Vf*h;
havg = .5*(hM+hM(mapP));
s1 = (-Drw*h + L*(havg.*nx))./J;
huM = Vf*hu;
huavg = .5*(huM+huM(mapP));
s2 = (-Drw*hu + L*(huavg.*nx))./J;
% uM = Vf*u;
% uavg = .5*(uM + uM(mapP));
% s2 = h.*(-Drw*u + L*(uavg.*nx))./J;

global a ep1 ep2
s1 = a.*ep1.*s1;
s2 = a.*ep2.*s2;
s1M = Vf*s1;
s2M = Vf*s2;
s1avg = .5*(s1M+s1M(mapP));
s2avg = .5*(s2M+s2M(mapP));
tau_ip = 1;
rhs1 = rhs1 + (-Drw*s1 + L*(s1avg.*nx + tau_ip*(hM(mapP)-hM)))./J;
rhs2 = rhs2 + (-Drw*s2 + L*(s2avg.*nx + tau_ip*(huM(mapP)-huM)))./J;

end
