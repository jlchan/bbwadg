clear
Globals2D

FinalTime = 1;

N = 6;
K1D = 10;
CFL = .5;
plotMesh = 0;

global tau
tau = .5;

%% set up physical domain 

Lx = 10; Ly = 5; ratiox = 1; ratioy = Ly/Lx;
[Nv, VX, VY, K, EToV] = unif_tri_mesh(round(ratiox*K1D),round(K1D*ratioy));
VX = (VX+1)*Lx; VY = VY*Ly;

StartUp2D;
BuildPeriodicMaps2D(max(VX)-min(VX),max(VY)-min(VY));
mapP = reshape(mapP,Nfp*Nfaces,K);

%% set up reference operators

global Vq Pq Lf Vff DNr DNs
global rxJ sxJ ryJ syJ nxJ nyJ 
global fxS1 fyS1 fxS2 fyS2 fxS3 fyS3 fxS4 fyS4
global pfun beta pavg plogmean vnormavg avg
global U1 U2 U3 U4

% plotting nodes
[rp sp] = EquiNodes2D(15); [rp sp] = xytors(rp,sp);

% volume quadrature nodes
[rq sq wq] = Cubature2D(2*N); 

% face quadrature nodes
[rq1D wq1D] = JacobiGQ(0,0,N);
rf = [rq1D; -rq1D; -ones(size(rq1D))];
sf = [-ones(size(rq1D)); rq1D; -rq1D];
wf = [wq1D; wq1D; wq1D];

% interpolation matrices
Vp = Vandermonde2D(N,rp,sp)/V;
Vq = Vandermonde2D(N,rq,sq)/V;
Vf = Vandermonde2D(N,rf,sf)/V;
Vq1D = Vandermonde1D(N,rq1D)/Vandermonde1D(N,JacobiGL(0,0,N));
Vff = kron(eye(3),Vq1D); % interpolates from face nodes to face quad points

% mass and projection matrices
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq)); % J's cancel out
Lf = M\(Vf'*diag(wf));

% scaled reference normals
z = zeros(size(rq1D));
e = ones(size(rq1D));
nrJ = [-z; e; -e];
nsJ = [-e; e; -z];

% "modal" differentiation matrices
Qr_modal = M*Dr;
Qs_modal = M*Ds;
Qr = Pq'*Qr_modal*Pq;
Qs = Pq'*Qs_modal*Pq;
E = Vf*Pq;
Br = diag(nrJ.*wf);
Bs = diag(nsJ.*wf);

% skew symmetric hybridized SBP operators
global QNr QNs
QNr = .5*[Qr-Qr' E'*Br;
    -Br*E zeros(length(wf))];
QNs = .5*[Qs-Qs' E'*Bs;
    -Bs*E zeros(length(wf))];

WN = diag([wq;wf]);
DNr = WN\QNr;
DNs = WN\QNs;

%% apply curved warping

a = 0/8;
x0 = Lx; y0 = 0;
x = x + Lx*a*cos(1/2*pi*(x-x0)/Lx).*cos(3/2*pi*(y-y0)/Ly);
y = y + Ly*a*sin(pi*(x-x0)/Lx).*cos(1/2*pi*(y-y0)/Ly);

%% make geometric terms

% get physical quadrature and plotting nodes
xq = Vq*x; yq = Vq*y;
xp = Vp*x; yp = Vp*y;

[rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);
rxJ = rx.*J;    sxJ = sx.*J;
ryJ = ry.*J;    syJ = sy.*J;

% compute scaled physical normals
nxJ = (Vf*rxJ).*nrJ + (Vf*sxJ).*nsJ;
nyJ = (Vf*ryJ).*nrJ + (Vf*syJ).*nsJ;
sJ = sqrt(nxJ.^2 + nyJ.^2);
nx = nxJ./sJ; ny = nyJ./sJ;

wJq = diag(wq)*(Vq*J);

if plotMesh
    rp1D = linspace(-1,1,100)';
    Vp1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,JacobiGL(0,0,N));
    Vfp = kron(eye(Nfaces),Vp1D);
    xfp = Vfp*x(Fmask(:),:);
    yfp = Vfp*y(Fmask(:),:);
    plot(xfp,yfp,'k.')
    hold on
    plot(x,y,'o')
    L2err = nan;
    return
end

%% fluxes
global gamma
gamma = 1.4;

rhoe = @(rho,rhou,rhov,E) E - .5*(rhou.^2+rhov.^2)./rho;
pcons = @(rho,rhou,rhov,E) (gamma-1)*rhoe(rho,rhou,rhov,E);
s = @(rho,rhou,rhov,E) log((gamma-1)*rhoe(rho,rhou,rhov,E)./(rho.^gamma));
V1 = @(rho,rhou,rhov,E) (-E + rhoe(rho,rhou,rhov,E).*(gamma + 1 - s(rho,rhou,rhov,E)))./(rhoe(rho,rhou,rhov,E));
V2 = @(rho,rhou,rhov,E) rhou./(rhoe(rho,rhou,rhov,E));
V3 = @(rho,rhou,rhov,E) rhov./(rhoe(rho,rhou,rhov,E));
V4 = @(rho,rhou,rhov,E) (-rho)./(rhoe(rho,rhou,rhov,E));

sV = @(V1,V2,V3,V4) gamma - V1 + (V2.^2+V3.^2)./(2*V4);
rhoeV  = @(V1,V2,V3,V4) ((gamma-1)./((-V4).^gamma)).^(1/(gamma-1)).*exp(-sV(V1,V2,V3,V4)/(gamma-1));
U1 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(-V4);
U2 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(V2);
U3 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(V3);
U4 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(1-(V2.^2+V3.^2)./(2*V4));

avg = @(x,y) .5*(x+y);
pfun = @(rho,u,v,E) (gamma-1)*(E-.5*rho.*(u.^2+v.^2));
beta = @(rho,u,v,E) rho./(2*pfun(rho,u,v,E));
pavg     = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER)     avg(rhoL,rhoR)./(2*avg(beta(rhoL,uL,vL,EL),beta(rhoR,uR,vR,ER)));
plogmean = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR)./(2*logmean(beta(rhoL,uL,vL,EL),beta(rhoR,uR,vR,ER)));
vnormavg = @(uL,vL,uR,vR) 2*(avg(uL,uR).^2 + avg(vL,vR).^2) - (avg(uL.^2,uR.^2) + avg(vL.^2,vR.^2));

% entropy potentials
psix = @(rho,rhou,rhov,E) (gamma-1)*rhou;
psiy = @(rho,rhou,rhov,E) (gamma-1)*rhov;

%% define initial condition using L2 projection

[rho u v p] = vortexSolution(x,y,0);
rhou = rho.*u;
rhov = rho.*v;
E    = p/(gamma-1) + .5*rho.*(u.^2+v.^2);

% Runge-Kutta residual storage
res1 = zeros(size(x));
res2 = zeros(size(x));
res3 = zeros(size(x));
res4 = zeros(size(x));

% compute time step size
CN = (N+1)*(N+2)/2; % based on trace constant
Jf = Vf*J;
CNh = max(CN*max(sJ(:)./Jf(:)));
dt = CFL*2/CNh;

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

rhstest = zeros(Nsteps,1);
entropy = zeros(Nsteps,1);

figure(1)
for i = 1:Nsteps
    for INTRK = 1:5
        
        rhoq  = Vq*rho;
        rhouq = Vq*rhou;
        rhovq = Vq*rhov;
        Eq    = Vq*E;
        
        % project to entropy variables for curvilinear        
        q1 = Pq*V1(rhoq,rhouq,rhovq,Eq);
        q2 = Pq*V2(rhoq,rhouq,rhovq,Eq);
        q3 = Pq*V3(rhoq,rhouq,rhovq,Eq);
        q4 = Pq*V4(rhoq,rhouq,rhovq,Eq);
                
        [rhs1 rhs2 rhs3 rhs4]  = RHS2Dsimple(q1,q2,q3,q4);        
        
        if (INTRK==5)
            r1 = (Vq*rhs1).*wJq.*(Vq*q1);
            r2 = (Vq*rhs2).*wJq.*(Vq*q2);
            r3 = (Vq*rhs3).*wJq.*(Vq*q3);
            r4 = (Vq*rhs4).*wJq.*(Vq*q4);
            rhstest(i) = rhstest(i) + sum(sum(r1+r2+r3+r4));
        end
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        res4 = rk4a(INTRK)*res4 + dt*rhs4;
        
        rho  = rho  + rk4b(INTRK)*res1;
        rhou = rhou + rk4b(INTRK)*res2;
        rhov = rhov + rk4b(INTRK)*res3;
        E    = E    + rk4b(INTRK)*res4;
        
    end;
    
    Sq = -rhoq.*s(rhoq,rhouq,rhovq,Eq);
    entropy(i) = sum(sum(wJq.*Sq));
    
    if  (mod(i,5)==0 || i==Nsteps)
        clf
        vv = Vp*rho;
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f, N = %d, K1D = %d',dt*i,N,K1D))
%         view(3)
        drawnow
    end    
end

[rq2 sq2 wq2] = Cubature2D(3*N);
Vq2 = Vandermonde2D(N,rq2,sq2)/V;
xq2 = Vq2*x; yq2 = Vq2*y;
wJq2 = diag(wq2)*(Vq2*J);
[rhoex uex vex pex] = vortexSolution(xq2,yq2,FinalTime);

rhouex = rhoex.*uex;
rhovex = rhoex.*vex;
Eex = pex/(gamma-1) + .5*rhoex.*(uex.^2+vex.^2);

rhoq = Vq2*rho;
rhouq = Vq2*rhou;
rhovq = Vq2*rhov;
Eq = Vq2*E;
err = wJq2.*((rhoq-rhoex).^2 + (rhouq-rhouex).^2 + (rhovq-rhovex).^2 + (Eq-Eex).^2);
L2err = sqrt(sum(err(:)));
fprintf('L2err = %8.8g\n',L2err)

dS = abs(entropy-entropy(1));
figure(2)
semilogy(dt*(1:Nsteps),dS,'o--');hold on
semilogy(dt*(1:Nsteps),abs(rhstest),'x--');hold on
legend('dS','rhstest')


function [rhs1 rhs2 rhs3 rhs4] = RHS2Dsimple(q1,q2,q3,q4)

Globals2D;

global Vq Pq Lf Vff DNr DNs
global rxJ sxJ ryJ syJ nxJ nyJ 
global fxS1 fyS1 fxS2 fyS2 fxS3 fyS3 fxS4 fyS4
global pfun beta pavg plogmean vnormavg avg
global U1 U2 U3 U4

%% compute jumps of entropy projected conservative variables

% evaluate at quad/surface points
q1q = Vq*q1; q2q = Vq*q2;  
q3q = Vq*q3; q4q = Vq*q4; 
rhoq  = U1(q1q,q2q,q3q,q4q); 
rhouq = U2(q1q,q2q,q3q,q4q); 
rhovq = U3(q1q,q2q,q3q,q4q); 
Eq    = U4(q1q,q2q,q3q,q4q); 
uq    = rhouq./rhoq; 
vq    = rhovq./rhoq; 
betaq = beta(rhoq,uq,vq,Eq);

% extract exterior entropy variable values 
q1M = q1(Fmask(:),:); 
q2M = q2(Fmask(:),:);
q3M = q3(Fmask(:),:); 
q4M = q4(Fmask(:),:);
q1P = q1M(mapP); 
q2P = q2M(mapP);
q3P = q3M(mapP); 
q4P = q4M(mapP);

% interpolate to quadrature points
q1M = Vff*q1M; q2M = Vff*q2M; q3M = Vff*q3M; q4M = Vff*q4M;
q1P = Vff*q1P; q2P = Vff*q2P; q3P = Vff*q3P; q4P = Vff*q4P;

rhoM   = U1(q1M,q2M,q3M,q4M);
rhouM  = U2(q1M,q2M,q3M,q4M);
rhovM  = U3(q1M,q2M,q3M,q4M);
EM     = U4(q1M,q2M,q3M,q4M);
uM     = rhouM./rhoM;
vM     = rhovM./rhoM;
betaM = beta(rhoM,uM,vM,EM);

rhoP  = U1(q1P,q2P,q3P,q4P);
rhouP = U2(q1P,q2P,q3P,q4P);
rhovP = U3(q1P,q2P,q3P,q4P);
EP    = U4(q1P,q2P,q3P,q4P);
uP    = rhouP./rhoP;
vP    = rhovP./rhoP;
betaP = beta(rhoP,uP,vP,EP);

% Lax-Friedrichs penalization
global gamma tau
unorm2 = (uM.^2+vM.^2);
p = (gamma-1)*(EM - .5*rhoM.*unorm2);
lamM = sqrt(unorm2)+sqrt(gamma*p./rhoM);
unorm2 = (uP.^2+vP.^2);
p = (gamma-1)*(EP - .5*rhoP.*unorm2);
lamP = sqrt(unorm2)+sqrt(gamma*p./rhoP);
LFc = max(lamP,lamM).*sJ;
Lf1 = tau*LFc.*(rhoP-rhoM);
Lf2 = tau*LFc.*(rhoP.*uP-rhoM.*uM);
Lf3 = tau*LFc.*(rhoP.*vP-rhoM.*vM);
Lf4 = tau*LFc.*(EP-EM);

[FxS1 FxS2 FxS3 FxS4 FyS1 FyS2 FyS3 FyS4] = euler_flux_2d(rhoM,rhoP,uM,uP,vM,vP,betaM,betaP);
fSf1 = nxJ.*FxS1 + nyJ.*FyS1;
fSf2 = nxJ.*FxS2 + nyJ.*FyS2;
fSf3 = nxJ.*FxS3 + nyJ.*FyS3;
fSf4 = nxJ.*FxS4 + nyJ.*FyS4;
fSf1 = fSf1  - .5*Lf1;
fSf2 = fSf2  - .5*Lf2;
fSf3 = fSf3  - .5*Lf3;
fSf4 = fSf4  - .5*Lf4;

%% compute flux differencing volume contributions

% stack vol/surface values
rhoN  = [rhoq; rhoM];
uN    = [uq; uM];
vN    = [vq; vM];
betaN = [betaq;betaM];

divF1 = zeros(size(DNr,1),K);
divF2 = zeros(size(DNr,1),K);
divF3 = zeros(size(DNr,1),K);
divF4 = zeros(size(DNr,1),K);
for e = 1:K
    
    [rhox, rhoy]   = meshgrid(rhoN(:,e));
    [ux, uy]       = meshgrid(uN(:,e));
    [vx, vy]       = meshgrid(vN(:,e));    
    [betax, betay] = meshgrid(betaN(:,e)); % no need to compute E 
    
    [FxS1 FxS2 FxS3 FxS4 FyS1 FyS2 FyS3 FyS4] = euler_flux_2d(rhox,rhoy,ux,uy,vx,vy,betax,betay);
        
    % assume constant geofacs - doesn't work properly if rxJ not constant!
    Qx = DNr.*rxJ(1,e) + DNs.*sxJ(1,e);
    Qy = DNr.*ryJ(1,e) + DNs.*syJ(1,e);       
   
    divF1(:,e) = sum(Qx.*FxS1,2) + sum(Qy.*FyS1,2);
    divF2(:,e) = sum(Qx.*FxS2,2) + sum(Qy.*FyS2,2);
    divF3(:,e) = sum(Qx.*FxS3,2) + sum(Qy.*FyS3,2);
    divF4(:,e) = sum(Qx.*FxS4,2) + sum(Qy.*FyS4,2);
end

rhs1 = 2*[Pq Lf]*divF1 + Lf*(fSf1);
rhs2 = 2*[Pq Lf]*divF2 + Lf*(fSf2);
rhs3 = 2*[Pq Lf]*divF3 + Lf*(fSf3);
rhs4 = 2*[Pq Lf]*divF4 + Lf*(fSf4);

rhs1 = -rhs1./J;
rhs2 = -rhs2./J;
rhs3 = -rhs3./J;
rhs4 = -rhs4./J;


end

function [rho u v p] = vortexSolution(x,y,t)

global gamma
x0 = 5;
y0 = 0;
beta = 5;
r2 = (x-x0-t).^2 + (y-y0).^2;

u = 1 - beta*exp(1-r2).*(y-y0)/(2*pi);
v = beta*exp(1-r2).*(x-x0-t)/(2*pi);
rho = 1 - (1/(8*gamma*pi^2))*(gamma-1)/2*(beta*exp(1-r2)).^2;
rho = rho.^(1/(gamma-1));
p = rho.^gamma;

if 0
    % pulse condition
    x0 = 5;
    rho = 2 + (abs(x-x0) < 2.5);
    u = 0*rho;
    v = 0*rho;
    p = rho.^gamma;    
end

end

function [FxS1 FxS2 FxS3 FxS4 FyS1 FyS2 FyS3 FyS4] = euler_flux_2d(rhoL,rhoR,uL,uR,vL,vR,betaL,betaR)

gamma = 1.4;

% optimized evaluations
rholog = logmean(rhoL,rhoR);
rhoavg = .5*(rhoL+rhoR);
uavg = .5*(uL+uR);
vavg = .5*(vL+vR);
vnavg = uL.*uR + vL.*vR; 
pa = rhoavg./(betaL+betaR);

EavgPavg = rholog./(2*(gamma-1)*logmean(betaL,betaR)) + pa + .5*rholog.*vnavg;

FxS1 = rholog.*uavg;      
FxS2 = FxS1.*uavg + pa;   
FxS3 = FxS1.*vavg;              
FxS4 = EavgPavg.*uavg;

FyS1 = rholog.*vavg;
FyS2 = FxS3;
FyS3 = FyS1.*vavg + pa;
FyS4 = EavgPavg.*vavg;

end