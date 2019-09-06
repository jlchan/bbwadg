clear

global N K

N = 7;
K = 32; % number of elements
CFL = .125; % should be O(1) - too large, soln blows up
FinalTime = 2;
tau = 1;

%% nodal DG

global Vf PN DNrskew Drw L mapP nx J

[r w] = JacobiGL(0,0,N);
rq = r; wq = w;
% [rq wq] = JacobiGQ(0,0,N);

V = Vandermonde1D(N,r);
Dr = GradVandermonde1D(N,r)/V;
Vq = Vandermonde1D(N,rq)/V;
Vf = Vandermonde1D(N,[-1;1])/V;

M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq));
L = M\(Vf');
Qr = diag(wq)*Vq*Dr*Pq;
Qrskew = (Qr-Qr');
QNr = [Qrskew Vf'*diag([-1;1])
    -diag([-1;1])*Vf zeros(2)];
WN = diag([wq;1;1]);
PN = M\[Vq' Vf']*WN;
Drw = M\(Qr');
DNrskew = WN\QNr;

VX = linspace(-1,1,K+1)';
h = VX(2)-VX(1);
x = zeros(N+1,K);
for e = 1:K
    x(:,e) = h*(1+r)/2 + VX(e);
end
J = h/2;

wJq = diag(wq)*repmat(J,N+1,K);
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
global fSn1fv fSn2fv
global fSn fS

% advection
fS = @(uL,uR) .5*(uL + uR);
fD = @(uL,uR) .5*(uR-uL);

% % Burgers
% fS = @(uL,uR) (uL.^2 + uL.*uR + uR.^2)/6;
% fD = @(uL,uR) .5*abs(.5*(uL+uR)).*(uR-uL);

% shallow water
% g = 9.81; % assume b = const for now
g = 1;
Sfun = @(h,u) .5*(h.*u.^2 + g*h.^2);
V1 = @(h,u) g*h - .5*u.^2;
V2 = @(h,u) u;
U1 = @(v1,v2) (v1+.5*v2.^2)/g;
U2 = @(v1,v2) v2;
fS1 = @(hL,hR,uL,uR) (hL.*uL + hR.*uR)/2;
fS2 = @(hL,hR,uL,uR) .5*(hL.*uL + hR.*uR).*.5.*(uL+uR) + .5*g.*(hL.*hR);
% fS2 = @(hL,hR,uL,uR) .5*(hL.*uL.^2 + hR.*uR.^2) + .5*g.*(hL.*hR);

% Lfc = @(hL,hR,uL,uR) abs(.5*(uL+uR))+sqrt(g*.5*(hL+hR)); % |u| + sqrt(gh)
Lf = @(h,u) abs(u)+sqrt(g*h);
Lfc = @(hL,hR,uL,uR) max(Lf(hL,uL),Lf(hR,uR)); % stronger LF
fD1 = @(hL,hR,uL,uR) .5*Lfc(hL,hR,uL,uR).*(hR-hL);
fD2 = @(hL,hR,uL,uR) .5*Lfc(hL,hR,uL,uR).*(hR.*uR-hL.*uL);

fSn1 = @(hL,hR,uL,uR,nx) fS1(hL,hR,uL,uR).*nx - tau*fD1(hL,hR,uL,uR);
fSn2 = @(hL,hR,uL,uR,nx) fS2(hL,hR,uL,uR).*nx - tau*fD2(hL,hR,uL,uR);

fSn1fv = @(hL,hR,uL,uR,nx) fS1(hL,hR,uL,uR).*nx - tau*fD1(hL,hR,uL,uR);
fSn2fv = @(hL,hR,uL,uR,nx) fS2(hL,hR,uL,uR).*nx - tau*fD2(hL,hR,uL,uR);

% interface flux
fSn = @(uL,uR,nx) fS(uL,uR).*nx - tau*fD(uL,uR); %.5*(uL.^b+uR.^b).*nx - tau*.5*(uR-uL).*abs(.5*(uR+uL)).^(b-1);

%% initial solution

h0 = @(x) rand(N+1,K)+1e-8; %
u0 = @(x) randn(N+1,K);

h0 = @(x) 1e-14 + exp(-50*x.^2); % gaussian initial condition
u0 = @(x) exp(-50*x.^2); % gaussian initial condition

% a = 1e-12;
% h0 = @(x) (1-a)*(x < 0) + a*ones(size(x));
u0 = @(x) zeros(size(x));

h = h0(x);
hu = h.*u0(x);

% S = h.^2 + h.*(hu./h).^2;

% plot(x,h,'o--');return

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
lam = abs(hu./h) + sqrt(g*h);
dt = CFL*min(dx(:)); %/max(lam(:)); % set time-step
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

figure(1)
resfv1 = zeros(N+1,K);
resfv2 = zeros(N+1,K);
res1 = zeros(N+1,K);
res2 = zeros(N+1,K);

global a ep1 ep2
for i = 1:Nsteps
    
    hh = h;
    hhuu = hu;
    hfv = h;
    hufv = hu;
    for INTRK = 1:3
        
        b = ssprkb(INTRK);
        
        a = ones(N+1,K);
        c = V\hh;
        sN = log(max(sum(c(N+1,:).^2,1)./sum(c.^2,1),sum(c(N,:).^2,1)./sum(c(1:N,:).^2,1)));
        ep1 = .5*J/N;
        ep2 = .5*J/N;
        s0 = log(1/N.^4);
        kappa = 10;
        sfun = .5*(1+sin(pi*(sN-s0)/(2*kappa)));
        a(:,sN < s0-kappa) = 0;
        a(:,(sN > s0-kappa) & (sN < s0 + kappa)) = repmat(sfun((sN > s0-kappa) & (sN < s0 + kappa)),N+1,1);
        a(:,(sN > s0+kappa)) = 1;
%         a(:) = 0;
        
        dry = min(hh) < 1e-13;
        
        hfvold = hfv;
        hufvold = hufv;
        [rhs1, rhs2] = rhs_fv(hfv,hufv);
        hfv = (1-b)*h + b*(hfv + dt*rhs1);
        hufv = (1-b)*hu + b*(hufv + dt*rhs2);
                
        [rhs1, rhs2] = rhs_DG(hh,hhuu);        
%         sum(sum(wJq.*(rhs1.*V1(hh,hhuu./hh) + rhs2.*V2(hh,hhuu./hh))))        
        hh = (1-b)*h + b*(hh + dt*rhs1);
        hhuu = (1-b)*hu + b*(hhuu + dt*rhs2);
        
        %         hfv = repmat(.5*wq'*hh,N+1,1);
        %         hufv = repmat(.5*wq'*hhuu,N+1,1);
        
        hh(isnan(hh(:))) = 0;
        hh(isinf(hh(:))) = 0;
        
        %         keyboard
        
        %         aa = max(min(-min(hfv)./(min(hh)-min(hfv)),1),0);
        aa = ones(1,K);
        aa(:,min(hh) < min(hfv) & min(hfv) < .05) = 0;
        aa(:,min(hh) < 1e-13) = 0;
        aa(:,max(abs(imag(hh)),[],1)>0) = 0;
        aa = repmat(aa,N+1,1);
%         aa(:) = 1;
        hh = hfv + aa.*(hh-hfv);
        hhuu = hufv + aa.*(hhuu-hufv);
        
        hh(:,dry) = 1e-14;
        hhuu(:,dry) = 0;
        
        if (min(hh(:))<1e-14 || min(hfv(:))<1e-14 || ~isreal(hfv))
            keyboard
        end
        
        
        
    end
    h = hh;
    hu = hhuu;
    
    
    if i==Nsteps || (mod(i,10)==0)
        %         plot(x,u,'bo--')
        plot(x,h,'bo-')
        hold on
        plot(x,hu,'r^-')
        plot(x,a,'ks')
        plot(x,aa,'m^')
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

function [rhsfv1 rhsfv2] = rhs_fv(h,hu)

global N K
global dx mapL mapR mapPfv fSn1 fSn2 nxfv

u = hu./h;

% 2nd order FV
dxL = .5*(dx + dx(mapL));
dxR = .5*(dx(mapR) + dx);

if 1
    hM = [h(:)';h(:)'];
    uM = [hu(:)';hu(:)']./hM;
else
    hsL = (h - h(mapL))./(dxL); % backwards diff slopes
    hsR = (h(mapR) - h)./(dxR); % forward diff slopes
    slopes = minmod([hsL(:)'; hsR(:)'])'; % minmod
    hM = [h(:) - slopes .* dx(:)/2, h(:) + slopes .* dx(:)/2]';
    
    husL = (hu - hu(mapL))./(dxL); % backwards diff slopes
    husR = (hu(mapR) - hu)./(dxR); % forward diff slopes
    slopes = minmod([husL(:)'; husR(:)'])'; % minmod
    uM = [hu(:) - slopes .* dx(:)/2, hu(:) + slopes .* dx(:)/2]'./hM;
end
hP = hM(mapPfv);
uP = uM(mapPfv);

% usL = (u - u(mapL))./(dxL); % backwards diff slopes
% usR = (u(mapR) - u)./(dxR); % forward diff slopes
% slopes = minmod([usL(:)'; usR(:)'])'; % minmod
% uL = u(:) - slopes .* dx(:)/2;
% uR = u(:) + slopes .* dx(:)/2;
% uM = [uL';uR'];
% uP = uM(mapPfv);


global fSn1fv fSn2fv
rhsfv1 = -reshape(sum(fSn1fv(hM,hP,uM,uP,nxfv),1),N+1,K)./dx;
rhsfv2 = -reshape(sum(fSn2fv(hM,hP,uM,uP,nxfv),1),N+1,K)./dx;


end

function [rhs1 rhs2] = rhs_DG(h,hu)

global N K
global Vf PN DNrskew Drw L mapP rxJ nx J
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
tau_ip = .5;
rhs1 = rhs1 + (-Drw*s1 + L*(s1avg.*nx + tau_ip*(hM(mapP)-hM)))./J;
rhs2 = rhs2 + (-Drw*s2 + L*(s2avg.*nx + tau_ip*(huM(mapP)-huM)))./J;

end
