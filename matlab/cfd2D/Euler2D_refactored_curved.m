clear
Globals2D

a = 0/8; % warping factor
FinalTime = 2;

N = 4;
K1D = 16;
CFL = .4;
plotMesh = 0;

global tau
tau = 1;

%% set up physical domain 

Lx = 7.5; Ly = 5; ratiox = 3/4; ratioy = .5;
% Lx = 10; Ly = 5; ratiox = 1; ratioy = Ly/Lx;
[Nv, VX, VY, K, EToV] = unif_tri_mesh(round(ratiox*K1D),round(K1D*ratioy));
VX = VX/max(abs(VX));  VY = VY/max(abs(VY));
VX = (VX+1)*Lx; VY = VY*Ly;

StartUp2D;
BuildPeriodicMaps2D(max(VX)-min(VX),max(VY)-min(VY));
mapP = reshape(mapP,Nfp*Nfaces,K);

%% set up reference operators

global M Vq Pq Lf Vff Vf Pfqf VqPq VqLf
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ wf

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

% reference normals
nrJ = [-zeros(size(rq1D)); ones(size(rq1D)); -ones(size(rq1D)); ];
nsJ = [-ones(size(rq1D)); ones(size(rq1D)); -zeros(size(rq1D)); ];

% "modal" differentiation matrices
Qr_modal = M*Dr;
Qs_modal = M*Ds;
Qr = Pq'*Qr_modal*Pq;
Qs = Pq'*Qs_modal*Pq;
E = Vf*Pq;
Br = diag(nrJ.*wf);
Bs = diag(nsJ.*wf);

VqLf = Vq*Lf;
VqPq = Vq*Pq;

% skew symmetric hybridized SBP operators
global DNr DNs WN 
QNr = .5*[Qr-Qr' E'*Br;
    -Br*E zeros(length(wf))];
QNs = .5*[Qs-Qs' E'*Bs;
    -Bs*E zeros(length(wf))];

WN = diag([wq;wf]);
DNr = WN\QNr;
DNs = WN\QNs;


%% make curvilinear mesh 

wadgProjEntropyVars = abs(a)>1e-8;

x0 = Lx; y0 = 0;
x = x + Lx*a*cos(1/2*pi*(x-x0)/Lx).*cos(3/2*pi*(y-y0)/Ly);
y = y + Ly*a*sin(pi*(x-x0)/Lx).*cos(1/2*pi*(y-y0)/Ly);

xq = Vq*x; yq = Vq*y;
xp = Vp*x; yp = Vp*y;

[rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Vq*Dr,Vq*Ds);
rxJ = rx.*J;    sxJ = sx.*J;
ryJ = ry.*J;    syJ = sy.*J;

[rx,sx,ry,sy,Jf] = GeometricFactors2D(x,y,Vf*Dr,Vf*Ds);
rxJf = rx.*Jf;    sxJf = sx.*Jf;
ryJf = ry.*Jf;    syJf = sy.*Jf;

nxJ = rxJf.*nrJ + sxJf.*nsJ;
nyJ = ryJf.*nrJ + syJf.*nsJ;
sJ = sqrt(nxJ.^2 + nyJ.^2);
nx = nxJ./sJ; ny = nyJ./sJ;

wJq = diag(wq)*(J);

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

global V1 V2 V3 V4
rhoe = @(rho,rhou,rhov,E) E - .5*(rhou.^2+rhov.^2)./rho;
pcons = @(rho,rhou,rhov,E) (gamma-1)*rhoe(rho,rhou,rhov,E);
s = @(rho,rhou,rhov,E) log((gamma-1)*rhoe(rho,rhou,rhov,E)./(rho.^gamma));
V1 = @(rho,rhou,rhov,E) (-E + rhoe(rho,rhou,rhov,E).*(gamma + 1 - s(rho,rhou,rhov,E)))./(rhoe(rho,rhou,rhov,E));
V2 = @(rho,rhou,rhov,E) rhou./(rhoe(rho,rhou,rhov,E));
V3 = @(rho,rhou,rhov,E) rhov./(rhoe(rho,rhou,rhov,E));
V4 = @(rho,rhou,rhov,E) (-rho)./(rhoe(rho,rhou,rhov,E));
% V1 = @(rho,rhou,rhov,E) (gamma-s(rho,rhou,rhov,E))/(gamma-1) - (rhou.^2+rhov.^2)./(rho.*2.*pcons(rho,rhou,rhov,E));
% V2 = @(rho,rhou,rhov,E) rhou./(pcons(rho,rhou,rhov,E));
% V3 = @(rho,rhou,rhov,E) rhov./(pcons(rho,rhou,rhov,E));
% V4 = @(rho,rhou,rhov,E) (-rho)./(pcons(rho,rhou,rhov,E));

global U1 U2 U3 U4
sV = @(V1,V2,V3,V4) gamma - V1 + (V2.^2+V3.^2)./(2*V4);
rhoeV  = @(V1,V2,V3,V4) ((gamma-1)./((-V4).^gamma)).^(1/(gamma-1)).*exp(-sV(V1,V2,V3,V4)/(gamma-1));
U1 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(-V4);
U2 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(V2);
U3 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(V3);
U4 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(1-(V2.^2+V3.^2)./(2*V4));

global fxS1 fyS1 fxS2 fyS2 fxS3 fyS3 fxS4 fyS4 psix psiy
global pfun beta pavg plogmean vnormavg avg

avg = @(x,y) .5*(x+y);
pfun = @(rho,u,v,E) (gamma-1)*(E-.5*rho.*(u.^2+v.^2));
beta = @(rho,u,v,E) rho./(2*pfun(rho,u,v,E));
pavg     = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER)     avg(rhoL,rhoR)./(2*avg(beta(rhoL,uL,vL,EL),beta(rhoR,uR,vR,ER)));
plogmean = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR)./(2*logmean(beta(rhoL,uL,vL,EL),beta(rhoR,uR,vR,ER)));
vnormavg = @(uL,vL,uR,vR) 2*(avg(uL,uR).^2 + avg(vL,vR).^2) - (avg(uL.^2,uR.^2) + avg(vL.^2,vR.^2));

fxS1 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(uL,uR);
fxS2 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(uL,uR).^2 + pavg(rhoL,uL,vL,EL,rhoR,uR,vR,ER);
fxS3 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(uL,uR).*avg(vL,vR);
fxS4 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) (plogmean(rhoL,uL,vL,EL,rhoR,uR,vR,ER)/(gamma-1) + pavg(rhoL,uL,vL,EL,rhoR,uR,vR,ER) + .5*logmean(rhoL,rhoR).*vnormavg(uL,vL,uR,vR)).*avg(uL,uR);

fyS1 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(vL,vR);
fyS2 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(uL,uR).*avg(vL,vR);
fyS3 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(vL,vR).^2 + pavg(rhoL,uL,vL,EL,rhoR,uR,vR,ER);
fyS4 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) (plogmean(rhoL,uL,vL,EL,rhoR,uR,vR,ER)/(gamma-1) + pavg(rhoL,uL,vL,EL,rhoR,uR,vR,ER) + .5*logmean(rhoL,rhoR).*vnormavg(uL,vL,uR,vR)).*avg(vL,vR);

% entropy potentials
psix = @(rho,rhou,rhov,E) (gamma-1)*rhou;
psiy = @(rho,rhou,rhov,E) (gamma-1)*rhov;

%% define initial condition using L2 projection

x0 = 0; y0 = 0;

[rhoq uq vq pq] = vortexSolution(xq,yq,0);
rho  = Pq*rhoq;
rhou = Pq*(rhoq.*uq);
rhov = Pq*(rhoq.*vq);
E    = Pq*(pq/(gamma-1) + .5*rhoq.*(uq.^2+vq.^2));

% Runge-Kutta residual storage
res1 = zeros(size(x));
res2 = zeros(size(x));
res3 = zeros(size(x));
res4 = zeros(size(x));

% compute time step size
CN = (N+1)*(N+2)/2; % guessing...
CNh = max(CN*max(sJ(:)./Jf(:)));
dt = CFL*2/CNh;

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

figure(1)
for i = 1:Nsteps
    for INTRK = 1:5
        
        rhoq  = Vq*rho;
        rhouq = Vq*rhou;
        rhovq = Vq*rhov;
        Eq    = Vq*E;
        
        % project to entropy variables for curvilinear
        if wadgProjEntropyVars
            q1 = Pq*((VqPq*(V1(rhoq,rhouq,rhovq,Eq).*J))./J);
            q2 = Pq*((VqPq*(V2(rhoq,rhouq,rhovq,Eq).*J))./J);
            q3 = Pq*((VqPq*(V3(rhoq,rhouq,rhovq,Eq).*J))./J);
            q4 = Pq*((VqPq*(V4(rhoq,rhouq,rhovq,Eq).*J))./J);
        else
            q1 = Pq*V1(rhoq,rhouq,rhovq,Eq);
            q2 = Pq*V2(rhoq,rhouq,rhovq,Eq);
            q3 = Pq*V3(rhoq,rhouq,rhovq,Eq);
            q4 = Pq*V4(rhoq,rhouq,rhovq,Eq);
        end
                
        [rhs1 rhs2 rhs3 rhs4]  = RHS2Dsimple(q1,q2,q3,q4);        
        
        if (INTRK==5)
            rhstest(i) = 0;
            if wadgProjEntropyVars                    
                for e = 1:K
                    r1 = M*((Vq'*diag(wq./J(:,e))*Vq)\(M*rhs1(:,e)));
                    r2 = M*((Vq'*diag(wq./J(:,e))*Vq)\(M*rhs2(:,e)));
                    r3 = M*((Vq'*diag(wq./J(:,e))*Vq)\(M*rhs3(:,e)));
                    r4 = M*((Vq'*diag(wq./J(:,e))*Vq)\(M*rhs4(:,e)));
                    rhstest(i) = rhstest(i) + sum((q1(:,e)'*(r1) + q2(:,e)'*(r2) + q3(:,e)'*(r3) + q4(:,e)'*(r4)));
                end
            else
                q1e = wJq.*(Vq*q1);
                q2e = wJq.*(Vq*q2);
                q3e = wJq.*(Vq*q3);
                q4e = wJq.*(Vq*q4);
                r1 = Vq*rhs1;
                r2 = Vq*rhs2;
                r3 = Vq*rhs3;
                r4 = Vq*rhs4;
                rhstest(i) = rhstest(i) + sum(sum(q1e.*r1 + q2e.*r2 + q3e.*r3 + q4e.*r4));
            end
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
        vv = real(Vp*rho);
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
wJq2 = diag(wq2)*(Vq2*Pq*J);
[rhoex uex vex pex] = vortexSolution(xq2,yq2,FinalTime);

rhouex = rhoex.*uex;
rhovex = rhoex.*vex;
% p = (gamma-1)*(E-.5*rho*(u^2+v^2));
Eex = pex/(gamma-1) + .5*rhoex.*(uex.^2+vex.^2);

rhoq = Vq2*rho;
rhouq = Vq2*rhou;
rhovq = Vq2*rhov;
Eq = Vq2*E;
err = wJq2.*((rhoq-rhoex).^2 + (rhouq-rhouex).^2 + (rhovq-rhovex).^2 + (Eq-Eex).^2);
L2err = sqrt(sum(err(:)));
fprintf('L2err = %8.8g\n',L2err)

dS = abs(entropy-entropy(1));
dS = entropy + max(abs(entropy));
figure(2)
semilogy(dt*(1:Nsteps),dS,'o--');hold on
semilogy(dt*(1:Nsteps),abs(rhstest),'x--');hold on
legend('dS','rhstest')


function [rhs1 rhs2 rhs3 rhs4] = RHS2Dsimple(q1,q2,q3,q4)

Globals2D;

global Vq Pq Lf Lff Vff Vf VqPq VqLf
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ nrJq nsJq
global mapPq
global fxS1 fyS1 fxS2 fyS2 fxS3 fyS3 fxS4 fyS4
global pfun beta pavg plogmean vnormavg avg
global U1 U2 U3 U4

global DNr DNs

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
betafq = beta(rhoM,uM,vM,EM);

rhoP  = U1(q1P,q2P,q3P,q4P);
rhouP = U2(q1P,q2P,q3P,q4P);
rhovP = U3(q1P,q2P,q3P,q4P);
EP    = U4(q1P,q2P,q3P,q4P);
uP    = rhouP./rhoP;
vP    = rhovP./rhoP;

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

fSf1 = nxJ.*fxS1(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + nyJ.*fyS1(rhoM,uM,vM,EM,rhoP,uP,vP,EP);
fSf2 = nxJ.*fxS2(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + nyJ.*fyS2(rhoM,uM,vM,EM,rhoP,uP,vP,EP);
fSf3 = nxJ.*fxS3(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + nyJ.*fyS3(rhoM,uM,vM,EM,rhoP,uP,vP,EP);
fSf4 = nxJ.*fxS4(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + nyJ.*fyS4(rhoM,uM,vM,EM,rhoP,uP,vP,EP);
fSf1 = fSf1  - .5*Lf1;
fSf2 = fSf2  - .5*Lf2;
fSf3 = fSf3  - .5*Lf3;
fSf4 = fSf4  - .5*Lf4;

%% compute flux differencing volume contributions

% stack vol/surface values
rhoN  = [rhoq; rhoM];
uN    = [uq; uM];
vN    = [vq; vM];
betaN = [betaq;betafq];

divF1 = zeros(size(DNr,1),K);
divF2 = zeros(size(DNr,1),K);
divF3 = zeros(size(DNr,1),K);
divF4 = zeros(size(DNr,1),K);
for e = 1:K
    
    [rhox, rhoy]   = meshgrid(rhoN(:,e));
    [ux, uy]       = meshgrid(uN(:,e));
    [vx, vy]       = meshgrid(vN(:,e));    
    [betax, betay] = meshgrid(betaN(:,e)); % no need to compute E 
    
    % optimized flux evaluations for compressible Euler
    rholog = logmean(rhox,rhoy);
    rhoavg = avg(rhox,rhoy);
    uavg = avg(ux,uy);
    vavg = avg(vx,vy);
    vnavg = 2*(uavg.^2 + vavg.^2) - (avg(ux.^2,uy.^2) + avg(vx.^2,vy.^2));
    pa = rhoavg./(2*avg(betax,betay));
    
    FxS1 = rholog.*uavg;      FyS1 = rholog.*vavg;
    FxS2 = FxS1.*uavg + pa;   FyS2 = FyS1.*uavg;
    FxS3 = FyS2;              FyS3 = FyS1.*vavg + pa;
    f4aux = rholog./(2*(gamma-1)*logmean(betax,betay)) + pa + .5*rholog.*vnavg;
    FxS4 = f4aux.*uavg;
    FyS4 = f4aux.*vavg;
    
    % split form for curved elements
    [rxJ1, rxJ2] = meshgrid([rxJ(:,e);rxJf(:,e)]);
    [sxJ1, sxJ2] = meshgrid([sxJ(:,e);sxJf(:,e)]);
    [ryJ1, ryJ2] = meshgrid([ryJ(:,e);ryJf(:,e)]);
    [syJ1, syJ2] = meshgrid([syJ(:,e);syJf(:,e)]);
    rxJK = avg(rxJ1,rxJ2);  sxJK = avg(sxJ1,sxJ2);
    ryJK = avg(ryJ1,ryJ2);  syJK = avg(syJ1,syJ2);
    
    Qx = DNr.*rxJK + DNs.*sxJK;
    Qy = DNr.*ryJK + DNs.*syJK;
    
    divF1(:,e) = sum(Qx.*FxS1,2) + sum(Qy.*FyS1,2);
    divF2(:,e) = sum(Qx.*FxS2,2) + sum(Qy.*FyS2,2);
    divF3(:,e) = sum(Qx.*FxS3,2) + sum(Qy.*FyS3,2);
    divF4(:,e) = sum(Qx.*FxS4,2) + sum(Qy.*FyS4,2);
end

rhs1 = 2*[VqPq VqLf]*divF1 + VqLf*(fSf1);
rhs2 = 2*[VqPq VqLf]*divF2 + VqLf*(fSf2);
rhs3 = 2*[VqPq VqLf]*divF3 + VqLf*(fSf3);
rhs4 = 2*[VqPq VqLf]*divF4 + VqLf*(fSf4);

% apply wadg for RHS 
rhs1 = -Pq*(rhs1./J);
rhs2 = -Pq*(rhs2./J);
rhs3 = -Pq*(rhs3./J);
rhs4 = -Pq*(rhs4./J);


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

% rho = (2 + sin(pi*(x - t)));
% u = ones(size(x));
% v = zeros(size(x));
% p = ones(size(x));

if 0
    % pulse condition
    x0 = 5;
    rho = 2 + (abs(x-x0) < 2.5);
    u = 0*rho;
    v = 0*rho;
    p = rho.^gamma;    
end

if 0 % sod
    x0 = 0;
    rhoL = 1; rhoR = .125;
    pL = 1; pR = .1;        
    
    rho = (rhoL*(x < x0) + rhoR*(x > x0));
    u = 0*x;
    v = 0*x;
    p = (pL*(x < x0) + pR*(x > x0));
    
end

end