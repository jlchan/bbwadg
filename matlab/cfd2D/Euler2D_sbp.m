% clear
Globals2D

N = 2;
K1D = 8;
FinalTime = 2;
CFL = .25;
global tau
tau = 1;
a = 0/16; % curv warping

Lx = 10; Ly = 5;
if 0    
    [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D,round(K1D*Ly/Lx));
    VX = VX/max(abs(VX));  VY = VY/max(abs(VY));
    VX = VX*Lx + Lx; VY = VY*Ly;
else
    [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
end

StartUp2D;
BuildPeriodicMaps2D(max(VX)-min(VX),max(VY)-min(VY));

% plotting nodes
[rp sp] = EquiNodes2D(15); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;
% plot(xp,yp,'o'); return

global Vq Pq Lq Lqf Vfqf Vfq Pfqf VqPq
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ nrJq nsJq
global mapPq
Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c

Vq = Vandermonde2D(N,rq,sq)/V;
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq)); % J's cancel out

VqPq = Vq*Pq;

xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

rxJ = rx.*J; sxJ = sx.*J;
ryJ = ry.*J; syJ = sy.*J;
rxJ = Vq*rxJ; sxJ = Vq*sxJ;
ryJ = Vq*ryJ; syJ = Vq*syJ;
J = Vq*J;

[rq1D wq1D] = JacobiGQ(0,0,ceil(Nq/2));
rfq = [rq1D; -rq1D; -ones(size(rq1D))];
sfq = [-ones(size(rq1D)); rq1D; -rq1D];
wfq = [wq1D; wq1D; wq1D];
Vq1D = Vandermonde1D(N,rq1D)/Vandermonde1D(N,JacobiGL(0,0,N));
% plot(rfq,sfq,'o')
Nfq = length(rq1D);

Nq = length(rq);

Vfq = Vandermonde2D(N,rfq,sfq)/V;
Vfqf = kron(eye(3),Vq1D);
Mf = Vfq'*diag(wfq)*Vfq;
Lq = M\(Vfq'*diag(wfq));

Pq1D = (Vq1D'*diag(wq1D)*Vq1D) \ (Vq1D'*diag(wq1D));
Pfqf = kron(eye(3),Pq1D);

nx = Vfqf*nx;
ny = Vfqf*ny;
sJ = Vfqf*sJ;
nxJ = (nx.*sJ);
nyJ = (ny.*sJ);
Fscale = Vfqf*Fscale;

nrJ = [-zeros(size(rq1D)); ones(size(rq1D)); -ones(size(rq1D)); ];
nsJ = [-ones(size(rq1D)); ones(size(rq1D)); -zeros(size(rq1D)); ];
Nq = length(rq);
nrJq = repmat(nrJ',Nq,1);
nsJq = repmat(nsJ',Nq,1);

% flux differencing operators
global Drq Dsq VfPq VqLq
Drq = (Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vfq*Pq);
Dsq = (Vq*Ds*Pq - .5*Vq*Lq*diag(nsJ)*Vfq*Pq);
VfPq = (Vfq*Pq);
VqLq = Vq*Lq;

DNr = [Drq .5*Vq*Lq*diag(nrJ);-.5*diag(nrJ)*Vfq*Pq .5*diag(nrJ)];
DNs = [Dsq .5*Vq*Lq*diag(nsJ);-.5*diag(nsJ)*Vfq*Pq .5*diag(nsJ)];
return

%% make quadrature face maps

xf = Vfq*x;
yf = Vfq*y;

mapMq = reshape(1:length(xf(:)),Nfq*Nfaces,K);
mapPq = mapMq;

for e = 1:K
    for f = 1:Nfaces
        enbr = EToE(e,f);
        if e ~= enbr % if it's a neighbor
            fnbr = EToF(e,f);
            id1 = (1:Nfq) + (f-1)*Nfq;
            id2 = (1:Nfq) + (fnbr-1)*Nfq;
            x1 = xf(id1,e); y1 = yf(id1,e);
            x2 = xf(id2,enbr); y2 = yf(id2,enbr);
            
            [X1 Y1] = meshgrid(x1,y1);
            [X2 Y2] = meshgrid(x2,y2);
            DX = (X1-X2').^2;
            DY = (Y1-Y2').^2;
            D = DX + DY;
            [p,~] = find(D<1e-8);
            
            if length(p) == 0
                % assume periodic boundary, find match in x,y
                [px,~] = find(DX<1e-8);
                [py,~] = find(DY<1e-8);
                if length(px)==0
                    p = py;
                elseif length(py)==0
                    p = px;
                else
                    keyboard
                end
                
            end
            mapPq(id1,e) = id2(p) + (enbr-1)*(Nfq*Nfaces);
        end
    end
end

%% make curvilinear mesh (still unstable?)

x0 = Lx; y0 = 0; 
x = x + Lx*a*cos(1/2*pi*(x-x0)/Lx).*cos(3/2*pi*(y-y0)/Ly);
y = y + Ly*a*sin(3/2*pi*(x-x0)/Lx).*cos(1/2*pi*(y-y0)/Ly);

xq = Vq*x; yq = Vq*y;
xp = Vp*x; yp = Vp*y;

if 0
    rp1D = linspace(-1,1,100)';
    Vp1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,JacobiGL(0,0,N));
    Vfp = kron(eye(Nfaces),Vp1D);
    xfp = Vfp*x(Fmask(:),:);
    yfp = Vfp*y(Fmask(:),:);
    plot(xfp,yfp,'k.')
    return
end

rxJ = zeros(Nq,K); sxJ = zeros(Nq,K);
ryJ = zeros(Nq,K); syJ = zeros(Nq,K);
J = zeros(Nq,K);
rxJf = zeros(Nfq*Nfaces,K); sxJf = zeros(Nfq*Nfaces,K);
ryJf = zeros(Nfq*Nfaces,K); syJf = zeros(Nfq*Nfaces,K);
Jf = zeros(Nfq*Nfaces,K);
for e = 1:K
    [rxk,sxk,ryk,syk,Jk] = GeometricFactors2D(x(:,e),y(:,e),Vq*Dr,Vq*Ds);
    rxJ(:,e) = rxk.*Jk;    sxJ(:,e) = sxk.*Jk;
    ryJ(:,e) = ryk.*Jk;    syJ(:,e) = syk.*Jk;
    J(:,e) = Jk;
    
    [rxk,sxk,ryk,syk,Jk] = GeometricFactors2D(x(:,e),y(:,e),Vfq*Dr,Vfq*Ds);
    rxJf(:,e) = rxk.*Jk;    sxJf(:,e) = sxk.*Jk;
    ryJf(:,e) = ryk.*Jk;    syJf(:,e) = syk.*Jk;
    Jf(:,e) = Jk;
end

nxJ = rxJf.*nrJ + sxJf.*nsJ;
nyJ = ryJf.*nrJ + syJf.*nsJ;

nx = nxJ./Jf;
ny = nyJ./Jf;
sJ = sqrt(nx.^2 + ny.^2);
nx = nx./sJ; ny = ny./sJ;
sJ = sJ.*Jf;

xf = Vfq*x;
yf = Vfq*y;
% plot(xf,yf,'o')
% hold on
% quiver(xf,yf,nxJ,nyJ)

%% fluxes
global gamma
gamma = 1.4;

rhoe = @(rho,rhou,rhov,E) E - .5*(rhou.^2+rhov.^2)./rho;
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


global fxS1 fyS1 fxS2 fyS2 fxS3 fyS3 fxS4 fyS4
global pfun beta pavg plogmean vnormavg avg
avg = @(x,y) .5*(x+y);
pfun = @(rho,u,v,E) (gamma-1)*(E-.5*rho.*(u.^2+v.^2));
beta = @(rho,u,v,E) rho./(2*pfun(rho,u,v,E));
pavg     = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER)     avg(rhoL,rhoR)./(2*avg(beta(rhoL,uL,vL,EL),beta(rhoR,uR,vR,ER)));
plogmean = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR)./(2*logmean(beta(rhoL,uL,vL,EL),beta(rhoR,uR,vR,ER)));
vnormavg = @(uL,vL,uR,vR) 2*(avg(uL,uR).^2 + avg(vL,vR).^2) - (avg(uL.^2,uR.^2) + avg(vL.^2,vR.^2));

% % shallow water test
% fxS4 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) avg(rhoL,rhoR).*avg(uL,uR);
% fxS2 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) avg(rhoL,rhoR).*avg(uL,uR).^2 + .5*avg(rhoL.^2,rhoR.^2);
% fxS3 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) avg(rhoL,rhoR).*avg(uL,uR).*avg(vL,vR);
% fxS1 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) 0*rhoL;
%
% fyS4 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) avg(rhoL,rhoR).*avg(vL,vR);
% fyS2 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) avg(rhoL,rhoR).*avg(uL,uR).*avg(vL,vR);
% fyS3 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) avg(rhoL,rhoR).*avg(vL,vR).^2 + .5*avg(rhoL.^2,rhoR.^2);
% fyS1 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) 0*rhoL;

fxS1 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(uL,uR);
fxS2 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(uL,uR).^2 + pavg(rhoL,uL,vL,EL,rhoR,uR,vR,ER);
fxS3 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(uL,uR).*avg(vL,vR);
fxS4 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) (plogmean(rhoL,uL,vL,EL,rhoR,uR,vR,ER)/(gamma-1) + pavg(rhoL,uL,vL,EL,rhoR,uR,vR,ER) + .5*logmean(rhoL,rhoR).*vnormavg(uL,vL,uR,vR)).*avg(uL,uR);

fyS1 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(vL,vR);
fyS2 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(uL,uR).*avg(vL,vR);
fyS3 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(vL,vR).^2 + pavg(rhoL,uL,vL,EL,rhoR,uR,vR,ER);
fyS4 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) (plogmean(rhoL,uL,vL,EL,rhoR,uR,vR,ER)/(gamma-1) + pavg(rhoL,uL,vL,EL,rhoR,uR,vR,ER) + .5*logmean(rhoL,rhoR).*vnormavg(uL,vL,uR,vR)).*avg(vL,vR);

% % conservative flux functions
% fxS1 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) rhoL.*uL;
% fxS2 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) rhoL.*uL.^2 + pfun(rhoL,uL,vL,EL);
% fxS3 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) rhoL.*uL.*vL;
% fxS4 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) (EL + (gamma-1)*(EL-.5*rhoL.*(uL.^2+vL.^2)./rhoL)).*rhoL.*uL;
%
% fyS1 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) rhoL.*vL;
% fyS2 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) rhoL.*uL.*vL;
% fyS3 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) rhoL.*vL.^2 + pfun(rhoL,uL,vL,EL);
% fyS4 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) (EL + pfun(rhoL,uL,vL,EL)./rhoL).*rhoL.*vL;


%% problem params setup

x0 = 0; y0 = 0;
[rhoq uq vq pq] = vortexSolution(xq,yq,0);
rho  = VqPq*rhoq;
rhou = VqPq*(rhoq.*uq);
rhov = VqPq*(rhoq.*vq);
E    = VqPq*(pq/(gamma-1) + .5*rhoq.*(uq.^2+vq.^2));


% vv = Vp*Pq*rho; color_line3(xp,yp,vv,vv,'.'); return

%%

global wJq
wJq = diag(wq)*(J);

% Runge-Kutta residual storage
res1 = zeros(Nq,K);
res2 = zeros(Nq,K);
res3 = zeros(Nq,K);
res4 = zeros(Nq,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
CNh = max(CN*max(sJ(:)./Jf(:)));
dt = CFL*2/CNh;

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

figure(1)
for i = 1:Nsteps
    for INTRK = 1:5
        
        rhoq  = rho;
        rhouq = rhou;
        rhovq = rhov;
        Eq    = E;
        
        % project to entropy variables
        q1 = Pq*V1(rhoq,rhouq,rhovq,Eq);
        q2 = Pq*V2(rhoq,rhouq,rhovq,Eq);
        q3 = Pq*V3(rhoq,rhouq,rhovq,Eq);
        q4 = Pq*V4(rhoq,rhouq,rhovq,Eq);
        
        % evaluate at quad/surface points
        q1q = Vq*q1;  q1M = Vfq*q1;
        q2q = Vq*q2;  q2M = Vfq*q2;
        q3q = Vq*q3;  q3M = Vfq*q3;
        q4q = Vq*q4;  q4M = Vfq*q4;
        rhoq  = U1(q1q,q2q,q3q,q4q); rhoM  = U1(q1M,q2M,q3M,q4M);
        rhouq = U2(q1q,q2q,q3q,q4q); rhouM = U2(q1M,q2M,q3M,q4M);
        rhovq = U3(q1q,q2q,q3q,q4q); rhovM = U3(q1M,q2M,q3M,q4M);
        Eq    = U4(q1q,q2q,q3q,q4q); EM    = U4(q1M,q2M,q3M,q4M);
        
        uq = rhouq./rhoq; uM = rhouM./rhoM;
        vq = rhovq./rhoq; vM = rhovM./rhoM;
        
        % extra LF flux info?
        QM{1} = rhoM;
        QM{2} = rhouM;  
        QM{3} = rhovM;
        QM{4} = EM;
        [rhs1 rhs2 rhs3 rhs4]  = RHS2D(rhoq,uq,vq,Eq,rhoM,uM,vM,EM,QM);
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        res4 = rk4a(INTRK)*res4 + dt*rhs4;
        
        rho  = rho  + rk4b(INTRK)*res1;
        rhou = rhou + rk4b(INTRK)*res2;
        rhov = rhov + rk4b(INTRK)*res3;
        E    = E    + rk4b(INTRK)*res4;
        
    end;
    
    Sq = -rho.*s(rho,rhou,rhov,E);
    S(i) = sum(sum(wJq.*Sq));
    
    if mod(i,10)==0 || i==Nsteps
        clf
        pp = Pq*rho;
        vv = real(Vp*pp);
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f',dt*i))
%         view(3)
        drawnow
    end
    
end

[rhoex uex vex pex] = vortexSolution(xq,yq,FinalTime);

err = wJq.*(rho-rhoex).^2;
L2err = sqrt(sum(err(:)))

figure(2)
dS = abs(S-S(1));
dS(end)
semilogy(dt*(1:Nsteps),dS(1:Nsteps),'--','linewidth',2)
hold on
ylabel('$\Delta U(t)$','fontsize',15,'Interpreter','latex')
xlabel('Time','fontsize',15)
set(gca,'fontsize',15)

function [rhs1 rhs2 rhs3 rhs4] = RHS2D(rhoq,uq,vq,Eq,rhoM,uM,vM,EM,QM)

Globals2D;

global Vq Pq Lq Lqf Vfqf Vfq Pfqf VqPq
% global Drq Dsq Lrq Lsq
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ nrJq nsJq
global mapPq
global fxS1 fyS1 fxS2 fyS2 fxS3 fyS3 fxS4 fyS4
global pfun beta pavg plogmean vnormavg avg
global Drq Dsq VfPq VqLq

Nq = size(Drq,2);

rhoP = rhoM(mapPq);
uP = uM(mapPq);
vP = vM(mapPq);
EP = EM(mapPq);

betaq = beta(rhoq,uq,vq,Eq);
betafq = beta(rhoM,uM,vM,EM);

% Lax-Friedrichs flux
global gamma
global tau
unorm2 = (uM.^2+vM.^2);
pM = (gamma-1)*(EM - .5*rhoM.*unorm2);
cvel = sqrt(gamma*pM./rhoM);
lam = sqrt(unorm2)+cvel;
LFc = max(lam(mapPq),lam);
QP{1} = QM{1}(mapPq); QP{2} = QM{2}(mapPq);
QP{3} = QM{3}(mapPq); QP{4} = QM{4}(mapPq);

% rhoUn = (rhoP.*uP-rhoM.*uM).*nx + (rhoP.*vP-rhoM.*vM).*ny;
dQ1 = QP{1}-QM{1};
dQ2 = QP{2}-QM{2};
dQ3 = QP{3}-QM{3};
dQ4 = QP{4}-QM{4};
Lf1 = tau*LFc.*dQ1.*sJ;
Lf2 = tau*LFc.*dQ2.*sJ;
Lf3 = tau*LFc.*dQ3.*sJ;
Lf4 = tau*LFc.*dQ4.*sJ;

fSf1 = nxJ.*fxS1(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + nyJ.*fyS1(rhoM,uM,vM,EM,rhoP,uP,vP,EP);
fSf2 = nxJ.*fxS2(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + nyJ.*fyS2(rhoM,uM,vM,EM,rhoP,uP,vP,EP);
fSf3 = nxJ.*fxS3(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + nyJ.*fyS3(rhoM,uM,vM,EM,rhoP,uP,vP,EP);
fSf4 = nxJ.*fxS4(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + nyJ.*fyS4(rhoM,uM,vM,EM,rhoP,uP,vP,EP);

fSf1 = fSf1  - .25*Lf1;
fSf2 = fSf2  - .25*Lf2;
fSf3 = fSf3  - .25*Lf3;
fSf4 = fSf4  - .25*Lf4;

rhs1 = zeros(Nq,K);
rhs2 = zeros(Nq,K);
rhs3 = zeros(Nq,K);
rhs4 = zeros(Nq,K);
for e = 1:K
    
    % local aritrhoMetic operations - form on the fly for GPU
    [rhox rhoy] = meshgrid(rhoq(:,e)); [rhofx rhofy] = meshgrid(rhoM(:,e),rhoq(:,e));
    [ux uy]     = meshgrid(uq(:,e));       [ufx ufy] = meshgrid(uM(:,e),uq(:,e));
    [vx vy]     = meshgrid(vq(:,e));       [vfx vfy] = meshgrid(vM(:,e),vq(:,e));
    [Ex Ey]     = meshgrid(Eq(:,e));       [Efx Efy] = meshgrid(EM(:,e),Eq(:,e));
    
    % avoiding geometric aliasing
    [rxJ1 rxJ2] = meshgrid(rxJ(:,e));  [sxJ1 sxJ2] = meshgrid(sxJ(:,e));
    [ryJ1 ryJ2] = meshgrid(ryJ(:,e));  [syJ1 syJ2] = meshgrid(syJ(:,e));
    rxJK = avg(rxJ1,rxJ2);  sxJK = avg(sxJ1,sxJ2);
    ryJK = avg(ryJ1,ryJ2);  syJK = avg(syJ1,syJ2);
    
    % optimized evaluations
    rholog = logmean(rhox,rhoy);
    rhoavg = avg(rhox,rhoy);
    uavg = avg(ux,uy);
    vavg = avg(vx,vy);
    %     vnavg = vnormavg(ux,vx,uy,vy); % 2*(avg(uL,uR).^2 + avg(vL,vR).^2) - (avg(uL.^2,uR.^2) + avg(vL.^2,vR.^2))
    vnavg = 2*(uavg.^2 + vavg.^2) - (avg(ux.^2,uy.^2) + avg(vx.^2,vy.^2));
    [betax betay] = meshgrid(betaq(:,e));
    pa = rhoavg./(2*avg(betax,betay));
    
    FxS1 = rholog.*uavg;      FyS1 = rholog.*vavg;
    FxS2 = FxS1.*uavg + pa;   FyS2 = FyS1.*uavg;
    FxS3 = FyS2;              FyS3 = FyS1.*vavg + pa;
    
    f4aux = rholog./(2*(gamma-1)*logmean(betax,betay)) + pa + .5*rholog.*vnavg;
    FxS4 = f4aux.*uavg;
    FyS4 = f4aux.*vavg;
        
    % premultiply by geofacs
    FrS1 = rxJK.*FxS1 + ryJK.*FyS1;     FsS1 = sxJK.*FxS1 + syJK.*FyS1;
    FrS2 = rxJK.*FxS2 + ryJK.*FyS2;     FsS2 = sxJK.*FxS2 + syJK.*FyS2;
    FrS3 = rxJK.*FxS3 + ryJK.*FyS3;     FsS3 = sxJK.*FxS3 + syJK.*FyS3;
    FrS4 = rxJK.*FxS4 + ryJK.*FyS4;     FsS4 = sxJK.*FxS4 + syJK.*FyS4;
    
    % bulk of GPU work: application of local operators
    divF1 = sum(Drq.*FrS1,2) + sum(Dsq.*FsS1,2);
    divF2 = sum(Drq.*FrS2,2) + sum(Dsq.*FsS2,2);
    divF3 = sum(Drq.*FrS3,2) + sum(Dsq.*FsS3,2);
    divF4 = sum(Drq.*FrS4,2) + sum(Dsq.*FsS4,2);
    
    %% flux/volume combos
        
    % optimized evaluations
    [betax betay] = meshgrid(betafq(:,e),betaq(:,e));
    rholog = logmean(rhofx,rhofy);
    rhoavg = avg(rhofx,rhofy);
    uavg = avg(ufx,ufy); vavg = avg(vfx,vfy);
    pa = rhoavg./(2*avg(betax,betay));
    vnavg = 2*(uavg.^2 + vavg.^2) - (avg(ufx.^2,ufy.^2) + avg(vfx.^2,vfy.^2));
    
    % averaging to prevent geometric aliasing
    [rxJ1 rxJ2] = meshgrid(rxJ(:,e),rxJf(:,e)); rxJK = avg(rxJ1,rxJ2)';
    f4aux = rholog./(2*(gamma-1)*logmean(betax,betay)) + pa + .5*rholog.*vnavg;
    
    FxSf1 = rholog.*uavg;       FySf1 = rholog.*vavg;
    FxSf2 = FxSf1.*uavg + pa;   FySf2 = FySf1.*uavg;
    FxSf3 = FySf2;              FySf3 = FySf1.*vavg + pa;
    
    FxSf4 = f4aux.*uavg;
    FySf4 = f4aux.*vavg;
    [sxJ1 sxJ2] = meshgrid(sxJ(:,e),sxJf(:,e)); sxJK = avg(sxJ1,sxJ2)';
    [ryJ1 ryJ2] = meshgrid(ryJ(:,e),ryJf(:,e)); ryJK = avg(ryJ1,ryJ2)';
    [syJ1 syJ2] = meshgrid(syJ(:,e),syJf(:,e)); syJK = avg(syJ1,syJ2)';
    
    Nq = size(FxSf1,1);
    FxSf1r = FxSf1.*nrJq; FxSf1s = FxSf1.*nsJq;
    FxSf2r = FxSf2.*nrJq; FxSf2s = FxSf2.*nsJq;
    FxSf3r = FxSf3.*nrJq; FxSf3s = FxSf3.*nsJq;
    FxSf4r = FxSf4.*nrJq; FxSf4s = FxSf4.*nsJq;
    
    FySf1r = FySf1.*nrJq; FySf1s = FySf1.*nsJq;
    FySf2r = FySf2.*nrJq; FySf2s = FySf2.*nsJq;
    FySf3r = FySf3.*nrJq; FySf3s = FySf3.*nsJq;
    FySf4r = FySf4.*nrJq; FySf4s = FySf4.*nsJq;
    
    FSf1 = rxJK.*FxSf1r + ryJK.*FySf1r + sxJK.*FxSf1s + syJK.*FySf1s;
    FSf2 = rxJK.*FxSf2r + ryJK.*FySf2r + sxJK.*FxSf2s + syJK.*FySf2s;
    FSf3 = rxJK.*FxSf3r + ryJK.*FySf3r + sxJK.*FxSf3s + syJK.*FySf3s;
    FSf4 = rxJK.*FxSf4r + ryJK.*FySf4r + sxJK.*FxSf4s + syJK.*FySf4s;
    
    if (e==round(K/2))
        keyboard
    end
    % bulk of GPU work: application of local operators
    fSproj1 = sum(VfPq.*FSf1',2);
    fSproj2 = sum(VfPq.*FSf2',2);
    fSproj3 = sum(VfPq.*FSf3',2);
    fSproj4 = sum(VfPq.*FSf4',2);
    
    % intermediate fluxes - form on the fly on GPU    
    f1 = FSf1 + repmat((fSf1(:,e) - fSproj1)',Nq,1) ;
    f2 = FSf2 + repmat((fSf2(:,e) - fSproj2)',Nq,1) ;
    f3 = FSf3 + repmat((fSf3(:,e) - fSproj3)',Nq,1) ;
    f4 = FSf4 + repmat((fSf4(:,e) - fSproj4)',Nq,1) ;
    
    %% project back to polynomial space - can incorporate WADG here
    
    rhs1(:,e) =  VqPq*(divF1 + .5*sum(VqLq.*f1,2));    
    rhs2(:,e) =  VqPq*(divF2 + .5*sum(VqLq.*f2,2));
    rhs3(:,e) =  VqPq*(divF3 + .5*sum(VqLq.*f3,2));
    rhs4(:,e) =  VqPq*(divF4 + .5*sum(VqLq.*f4,2));
    
end

rhs1 = -2*VqPq*(rhs1./J);
rhs2 = -2*VqPq*(rhs2./J);
rhs3 = -2*VqPq*(rhs3./J);
rhs4 = -2*VqPq*(rhs4./J);

end



function [rho u v p] = vortexSolution(x,y,t)

global gamma
x0 = 5;
y0 = 0;
beta = 5;
r = sqrt((x-x0-t).^2 + (y-y0).^2);

u = 1 - beta*exp(1-r.^2).*(y-y0)/(2*pi);
v = beta*exp(1-r.^2).*(x-x0-t)/(2*pi);
% rho = 1 - (gamma-1)/(16*gamma*pi^2)*beta^2*exp(2*(1-r.^2));
rho = 1 - (1/(8*gamma*pi^2))*(gamma-1)/2*(beta*exp(1-r.^2)).^2;
rho = rho.^(1/(gamma-1));
p = rho.^gamma;

% rho = (2 + sin(pi*(x - t)));
% u = ones(size(x));
% v = zeros(size(x));
% p = ones(size(x));

% rho = ones(size(x));
% u = ones(size(x));
% v = zeros(size(x));
% p = ones(size(x));

% pulse condition
x0 = 0; y0 = 0;
rho = 2 + (abs(x-x0) < .5).*(abs(y-y0) < .5);
u = 0*rho;
v = u;
p = rho.^gamma;

end
