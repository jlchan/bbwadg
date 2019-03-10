clear
clear -globals
Globals2D

% aVX = .125; % vertex warping factor
a = 0*.125; % warping factor

N = 3;
K1D = 20;

Mgeo = N;
[rq1D_face wq1D_face] = JacobiGQ(0,0,N);

FinalTime = 6.0;

wadgProjEntropyVars = abs(a)>1e-8;
CFL = .75;
global tau
tau = 1;
plotMesh = 0;

%%

Lx = 1; Ly = 3/4; ratiox = 1; ratioy = Ly/Lx;
[Nv, VX, VY, K, EToV] = QuadMesh2D(round(ratiox*K1D),round(ratioy*K1D));
% VX = VX + aVX*cos(pi/2*VX).*sin(pi*VY);
% VY = VY + aVX*sin(pi*VX).*cos(pi/2*VY);

% ids = abs(abs(VX)-1)>1e-8 & abs(abs(VY)-1)>1e-8;
% VY(ids) = VY(ids) + .05*randn(size(VX(ids)));
% VX = VX/max(abs(VX));  VY = VY/max(abs(VY));
VX = VX*Lx; VY = VY*Ly;

% VX = VX + .75*randn(size(VX));
% VY = VY + .75*randn(size(VX));
StartUp2D;
BuildPeriodicMaps2D(max(VX)-min(VX),max(VY)-min(VY));
% BuildPeriodicMaps2D(max(VX)-min(VX),0);

global M Vq Pq Lq Vfqf Vfq Pfqf VqPq VqLq
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ wfq
global mapPq

% vol nodes
[rq1D wq1D] = JacobiGQ(0,0,N);
[rq sq] = meshgrid(rq1D);
rq = rq(:); sq = sq(:);
[wr ws] = meshgrid(wq1D);
wq = wr(:).*ws(:);

Vq = Vandermonde2D(N,rq,sq)/V;
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq)); % J's cancel out

% face nodes
rq1D = rq1D_face;
wq1D = wq1D_face;
%[rq1D wq1D] = JacobiGQ(0,0,ceil((N-1)/2)); % min for entropy stability on affine meshes
% rq1D = rq1D*(1-1e-10);
% rq1D = [-1+(1+rq1D)/2; (1+rq1D)/2];
% wq1D = [wq1D; wq1D]/2;
% rq1D = rq1D*(1-5e-10);
Nfq = length(rq1D);

e = ones(size(rq1D));
rfq = [rq1D; e; rq1D; -e];
sfq = [-e; rq1D; e; rq1D];
wfq = [wq1D; wq1D; wq1D; wq1D];
V1D = Vandermonde1D(N,JacobiGL(0,0,N));
Vq1D = Vandermonde1D(N,rq1D)/V1D;
Vfqf = kron(eye(Nfaces),Vq1D);
Vfq = Vandermonde2D(N,rfq,sfq)/V;
Mf = Vfq'*diag(wfq)*Vfq;
Lq = M\(Vfq'*diag(wfq));

nrJ = [0*e; e; 0*e; -e]; % sJ = 2 for all faces,
nsJ = [-e; 0*e; e; 0*e];

% flux differencing operators
VfPq = (Vfq*Pq);
VqLq = Vq*Lq;
VqPq = Vq*Pq;

global DNr DNs WN
Drq = (Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vfq*Pq);
Dsq = (Vq*Ds*Pq - .5*Vq*Lq*diag(nsJ)*Vfq*Pq);
DNr = [Drq .5*VqLq*diag(nrJ);
    -.5*diag(nrJ)*VfPq .5*diag(nrJ)];
DNs = [Dsq .5*VqLq*diag(nsJ);
    -.5*diag(nsJ)*VfPq .5*diag(nsJ)];
WN = diag([wq;wfq]);

% weak derivative operators
Drqw = (Vq*(M\(Dr'*Vq'*diag(wq))) - .5*Vq*Lq*diag(nrJ)*Vfq*Pq);
Dsqw = (Vq*(M\(Ds'*Vq'*diag(wq))) - .5*Vq*Lq*diag(nsJ)*Vfq*Pq);

DNrw = [Drqw -.5*VqLq*diag(nrJ);
    .5*diag(nrJ)*VfPq .5*diag(nrJ)];
DNsw = [Dsqw -.5*VqLq*diag(nsJ);
    .5*diag(nsJ)*VfPq .5*diag(nsJ)];

QNrskew = .5*(WN*DNr - (WN*DNr)');
QNsskew = .5*(WN*DNs - (WN*DNs)');

% make skew symmetric diff matrices
DNr = diag(1./[wq;wfq])*QNrskew;
DNs = diag(1./[wq;wfq])*QNsskew;

DNr(abs(DNr)<1e-8) = 0;
DNs(abs(DNs)<1e-8) = 0;

DNr = sparse(DNr);
DNs = sparse(DNs);


QNr = WN*[Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vfq*Pq .5*VqLq*diag(nrJ);
    -.5*diag(nrJ)*VfPq .5*diag(nrJ)];
BNr = [0*Drq 0*VqLq;
    0*VfPq diag(wfq.*nrJ)];

rN = [rq;rfq];
sN = [sq;sfq];
u = [Vq;Vfq]*randn(Np,1);
v = sN.^(N-1)+rN.^(N-1);

norm(BNr-(QNr+QNr'))
norm(v'*QNr*u - v'*(BNr-QNr')*u)

% plotting nodes
[rp sp] = EquiNodes2D(15); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

%% make quadrature face maps

xf = Vfq*x;
yf = Vfq*y;

mapMq = reshape(1:length(xf(:)),Nfq*Nfaces,K);
mapPq = mapMq;

tol = 1e-10;
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
            [p,~] = find(D<tol);
            
            % NOTE - does not work if K1D is too small!!
            if length(p) == 0
                % assume periodic boundary, find match in x,y
                [px,~] = find(DX<tol);
                [py,~] = find(DY<tol);
                if length(px)==0
                    p = py;
                elseif length(py)==0
                    p = px;
                else
                    keyboard
                end
            end
            fids = id2(p) + (enbr-1)*(Nfq*Nfaces);
            mapPq(id1,e) = fids(:);
        end
    end
end

global mapBwall
mapBwall = find(abs(yf-max(y(:)))<1e-8 | abs(yf-min(y(:)))<1e-8);


%% make curvilinear mesh

x0 = Lx; y0 = 0;
% x0 = 0; y0 = 0; Lx = 1; Ly = 1;
dx = cos(1/2*pi*(x-x0)/Lx).*cos(3/2*pi*(y-y0)/Ly);
x = x + Lx*a*dx;
dy = sin(4/2*pi*(x-x0)/Lx).*cos(1/2*pi*(y-y0)/Ly);
y = y + Ly*a*dy;

% filter nodes by interpolating down to degree M and
[rM sM] = meshgrid(JacobiGL(0,0,Mgeo)); rM = rM(:); sM = sM(:);
F = (Vandermonde2DQuad(Mgeo,r,s)/Vandermonde2DQuad(Mgeo,rM,sM))*(Vandermonde2DQuad(N,rM,sM)/Vandermonde2DQuad(N,r,s));
x = F*x;
y = F*y;

xq = Vq*x; yq = Vq*y;
xp = Vp*x; yp = Vp*y;
xf = Vfq*x;    yf = Vfq*y;

[rx sx ry sy J] = GeometricFactors2D(x,y,Vq*Dr,Vq*Ds);
rxJ = rx.*J; sxJ = sx.*J;
ryJ = ry.*J; syJ = sy.*J;

% rxJf = Vfq*rxJ; sxJf = Vfq*sxJ;
% ryJf = Vfq*ryJ; syJf = Vfq*syJ;

[rxf,sxf,ryf,syf,Jf] = GeometricFactors2D(x,y,Vfq*Dr,Vfq*Ds);
rxJf = rxf.*Jf; sxJf = sxf.*Jf;
ryJf = ryf.*Jf; syJf = syf.*Jf;

nxJ = rxJf.*nrJ + sxJf.*nsJ;
nyJ = ryJf.*nrJ + syJf.*nsJ;

% nx = nxJ./Jf;
% ny = nyJ./Jf;
sJ = sqrt(nxJ.^2 + nyJ.^2);
% nx = nxJ./sJ; ny = nyJ./sJ;
% sJ = sJ.*Jf;

if plotMesh
    rp1D = linspace(-1,1,100)';
    Vp1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,JacobiGL(0,0,N));
    Vfp = kron(eye(Nfaces),Vp1D);
    xfp = Vfp*x(Fmask(:),:);
    yfp = Vfp*y(Fmask(:),:);
    plot(xfp,yfp,'k.')
    hold on
    %plot(x,y,'b.','markersize',14)
    for e = 1:K
        xe = reshape(x(:,e),N+1,N+1);
        ye = reshape(y(:,e),N+1,N+1);
        for i = 1:N+1
            plot(Vp1D*xe(:,i),Vp1D*ye(:,i),'k-');
            plot(Vp1D*(xe(i,:)'),Vp1D*(ye(i,:)'),'k-');
        end
    end
    axis equal
    L2err = nan;
    axis off
    return
end

%% fluxes
global gamma
gamma = 1.4;

global V1 V2 V3 V4
rhoe = @(rho,rhou,rhov,E) E - .5*(rhou.^2+rhov.^2)./rho;
pcons = @(rho,rhou,rhov,E) (gamma-1)*rhoe(rho,rhou,rhov,E);
sfun = @(rho,rhou,rhov,E) log((gamma-1)*rhoe(rho,rhou,rhov,E)./(rho.^gamma));
V1 = @(rho,rhou,rhov,E) (-E + rhoe(rho,rhou,rhov,E).*(gamma + 1 - sfun(rho,rhou,rhov,E)))./(rhoe(rho,rhou,rhov,E));
V2 = @(rho,rhou,rhov,E) rhou./(rhoe(rho,rhou,rhov,E));
V3 = @(rho,rhou,rhov,E) rhov./(rhoe(rho,rhou,rhov,E));
V4 = @(rho,rhou,rhov,E) (-rho)./(rhoe(rho,rhou,rhov,E));

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


%% problem params setup

x0 = 0; y0 = 0;

[rhoq uq vq pq] = vortexSolution(xq,yq,0);
rho  = Vq*Pq*rhoq;
rhou = Vq*Pq*(rhoq.*uq);
rhov = Vq*Pq*(rhoq.*vq);
E    = Vq*Pq*(pq/(gamma-1) + .5*rhoq.*(uq.^2+vq.^2));

% vv = (rhou.^2+rhov.^2)./rho.^2;
% color_line3(xq,yq,vv,vv,'.')
% return

%%

global wJq
wJq = diag(wq)*(J);

global VqPN
VqPN = 2*[VqPq VqLq];

% Runge-Kutta residual storage
Nq = length(rq);
res1 = zeros(Nq,K);
res2 = zeros(Nq,K);
res3 = zeros(Nq,K);
res4 = zeros(Nq,K);

% compute time step size
CN = (N+1)*(N+2)/2; % for quads
CNh = max(CN*max(sJ(:)./Jf(:)));
dt = CFL*1/CNh;

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

q1 = V1(rho,rhou,rhov,E);
q2 = V2(rho,rhou,rhov,E);
q3 = V3(rho,rhou,rhov,E);
q4 = V4(rho,rhou,rhov,E);

interval = 5;
Usnap = zeros(4*Np*K,ceil(Nsteps/interval));
Vsnap = zeros(4*Np*K,ceil(Nsteps/interval));
sk = 1;

figure(1)
for i = 1:Nsteps
    for INTRK = 1:5
        
        q1q = V1(rho,rhou,rhov,E);
        q2q = V2(rho,rhou,rhov,E);
        q3q = V3(rho,rhou,rhov,E);
        q4q = V4(rho,rhou,rhov,E);
        
        q1M = VfPq*q1q;
        q2M = VfPq*q2q;
        q3M = VfPq*q3q;
        q4M = VfPq*q4q;
        
        rhoq  = U1(q1q,q2q,q3q,q4q); rhoM  = U1(q1M,q2M,q3M,q4M);
        rhouq = U2(q1q,q2q,q3q,q4q); rhouM = U2(q1M,q2M,q3M,q4M);
        rhovq = U3(q1q,q2q,q3q,q4q); rhovM = U3(q1M,q2M,q3M,q4M);
        Eq    = U4(q1q,q2q,q3q,q4q); EM    = U4(q1M,q2M,q3M,q4M);
        
        uq = rhouq./rhoq; uM = rhouM./rhoM;
        vq = rhovq./rhoq; vM = rhovM./rhoM;
        
        % extra LF flux info
        [rhs1 rhs2 rhs3 rhs4]  = RHS2Dsimple(rhoq,uq,vq,Eq,rhoM,uM,vM,EM);
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        res4 = rk4a(INTRK)*res4 + dt*rhs4;
        
        rho  = rho  + rk4b(INTRK)*res1;
        rhou = rhou + rk4b(INTRK)*res2;
        rhov = rhov + rk4b(INTRK)*res3;
        E    = E    + rk4b(INTRK)*res4;
        
    end;
    
    if mod(i,interval)==0 % t > 1.001652892561984
        q1 = V1(rho,rhou,rhov,E);
        q2 = V2(rho,rhou,rhov,E);
        q3 = V3(rho,rhou,rhov,E);
        q4 = V4(rho,rhou,rhov,E);
        Usnap(:,sk) = [rho(:);rhou(:);rhov(:);E(:)];
        Vsnap(:,sk) = [q1(:);q2(:);q3(:);q4(:)];
        tsnap(sk) = i*dt;
        sk = sk + 1;
    end
    
    %     Sq = -rho.*s(rho,rhou,rhov,E);
    %     entropy(i) = sum(sum(wJq.*Sq));
    
    if  (mod(i,5)==0 || i==Nsteps)
        clf
        vv = Vp*Pq*(rhou); %Vp*Pq*((rhou./rho).^2+(rhov./rho).^2);
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f, N = %d, K1D = %d, tstep %d/%d',dt*i,N,K1D,i,Nsteps))
        drawnow
    end
end

return

%%

Ntotal = size(Usnap,1)/4;
ids = 1:Ntotal;

Us = [Usnap(ids,:) Usnap(ids+Ntotal,:) Usnap(ids+2*Ntotal,:) Usnap(ids+3*Ntotal,:)];
Vs = [Vsnap(ids,:) Vsnap(ids+Ntotal,:) Vsnap(ids+2*Ntotal,:) Vsnap(ids+3*Ntotal,:)];
[U, S, ~] = svd([Us Vs],0);
[Uh,Sh,~] = svd([U(:,1:Nmodes) ones(Ntotal,1)],0);

%%
Nmodes = 25;

for i = 1:5:size(Usnap,2)-1    
    U1 = reshape(Usnap(:,i),Np*K,4);
    rho = reshape(U1(:,1),Np,K);
    rhou = reshape(U1(:,2),Np,K);
    rhov = reshape(U1(:,3),Np,K);
    E = reshape(U1(:,4),Np,K);
    
%     vv = rhou;
    vv = reshape(Uh*(Uh'*rhou(:)),Np,K);
    
    %     figure
    clf    
    vv = Vp*Pq*vv;
    color_line3(xp,yp,vv,vv,'.');
    axis equal
    axis tight
    colorbar
    drawnow
end

%% compute deim quadrature

Nmodes = 100;
Nq = Nmodes;
UDEIM(:,1) = Uh(:,1);
[~,id] = max(abs(Uh(:,1)));
p = id;
for j = 2:Nq
    r = Uh(:,j)-UDEIM*(UDEIM(p,:)\Uh(p,j));
    [~,id] = max(abs(r));
    p(j) = id;
    UDEIM = [UDEIM r];    
end
xDEIM = xq(p);
yDEIM = yq(p);
wDEIM = wJq(:)'*(UDEIM/UDEIM(p,:));
% wDEIM2 = Uh(p,1:Nq)'\(Uh(:,1:Nq)'*wJq(:));

plot(xq,yq,'b.')
hold on
plot(xDEIM,yDEIM,'ro','markersize',16,'linewidth',2)

clf
hold on
vv = UDEIM/UDEIM(p,:);
vv = vv(:,1);
color_line3(xq(:),yq(:),vv,vv,'.')
plot3(xDEIM(:),yDEIM(:),vv(p),'ro','markersize',16,'linewidth',2)

%%

function [rhs1 rhs2 rhs3 rhs4] = RHS2Dsimple(rhoq,uq,vq,Eq,rhoM,uM,vM,EM,QM)

Globals2D;

global M Vq Pq Lq Lqf Vfqf Vfq VqPq VqLq
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ
global mapPq
global fxS1 fyS1 fxS2 fyS2 fxS3 fyS3 fxS4 fyS4
global pfun beta pavg plogmean vnormavg avg
global Drq Dsq VfPq

global DNr DNs

rhoP = rhoM(mapPq);
uP = uM(mapPq);
vP = vM(mapPq);
EP = EM(mapPq);

% walls = top/bottom
global mapBwall
uP(mapBwall) = uM(mapBwall);
vP(mapBwall) = -vM(mapBwall);
rhoP(mapBwall) = rhoM(mapBwall);
EP(mapBwall) = EM(mapBwall);

betaq = beta(rhoq,uq,vq,Eq);
betafq = beta(rhoM,uM,vM,EM);

% Lax-Friedrichs flux
global gamma tau
unorm2 = (uM.^2+vM.^2);
pM = (gamma-1)*(EM - .5*rhoM.*unorm2);
cvel = sqrt(gamma*pM./rhoM);
lam = sqrt(unorm2)+cvel;
LFc = max(lam(mapPq),lam);

dQ1 = rhoP-rhoM;
dQ2 = rhoP.*uP-rhoM.*uM;
dQ3 = rhoP.*vP-rhoM.*vM;
dQ4 = EP-EM;
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

rho = [rhoq; rhoM];
u = [uq; uM];
v = [vq; vM];
% E = [Eq; EM];
betaN = [betaq;betafq];

divF1 = zeros(size(DNr,1),K);
divF2 = zeros(size(DNr,1),K);
divF3 = zeros(size(DNr,1),K);
divF4 = zeros(size(DNr,1),K);
for e = 1:K
    
    [rhox, rhoy] = meshgrid(rho(:,e));
    [ux, uy] = meshgrid(u(:,e));
    [vx, vy] = meshgrid(v(:,e));
    %     [Ex, Ey] = meshgrid(E(:,e)); % no need - used in beta
    [betax, betay] = meshgrid(betaN(:,e));
    
    % optimized evaluations
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
    
    %     [rxJ1, rxJ2] = meshgrid([rxJ(:,e);rxJf(:,e)]);
    %     [sxJ1, sxJ2] = meshgrid([sxJ(:,e);sxJf(:,e)]);
    %     [ryJ1, ryJ2] = meshgrid([ryJ(:,e);ryJf(:,e)]);
    %     [syJ1, syJ2] = meshgrid([syJ(:,e);syJf(:,e)]);
    %     rxJK = avg(rxJ1,rxJ2);  sxJK = avg(sxJ1,sxJ2);
    %     ryJK = avg(ryJ1,ryJ2);  syJK = avg(syJ1,syJ2);
    rxJK = rxJ(1,e); ryJK = ryJ(1,e);
    sxJK = sxJ(1,e); syJK = syJ(1,e);
    
    Dx = DNr.*rxJK + DNs.*sxJK;
    Dy = DNr.*ryJK + DNs.*syJK;
    
    divF1(:,e) = sum(Dx.*FxS1,2) + sum(Dy.*FyS1,2);
    divF2(:,e) = sum(Dx.*FxS2,2) + sum(Dy.*FyS2,2);
    divF3(:,e) = sum(Dx.*FxS3,2) + sum(Dy.*FyS3,2);
    divF4(:,e) = sum(Dx.*FxS4,2) + sum(Dy.*FyS4,2);
end

global VqPN
rhs1 = VqPN*divF1 + VqLq*(fSf1);
rhs2 = VqPN*divF2 + VqLq*(fSf2);
rhs3 = VqPN*divF3 + VqLq*(fSf3);
rhs4 = VqPN*divF4 + VqLq*(fSf4);

% apply wadg at the end
rhs1 = -(rhs1./J);
rhs2 = -(rhs2./J);
rhs3 = -(rhs3./J);
rhs4 = -(rhs4./J);


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


opt = 1;
if opt==1 % jet
    rho = ones(size(x));
    p = ones(size(x));
    
    width = .1;
    sig = 100;
    ff = 1-(1./(1+exp(-sig*(y-width))) + 1./(1+exp(sig*(y+width))));
    dv = .1*sin(2*pi*x).*ff;
    du = .5*ff;
    v = zeros(size(x)) + dv;
    u = -.5*ones(size(x)) + du;
    
elseif opt==2
    % pulse condition
    x0 = 5;
    rho = 2 + (abs(x-x0) < 2.5);
    u = 0*rho;
    v = 0*rho;
    p = rho.^gamma;
elseif opt==3 %sod
    x0 = 0;
    rhoL = 1; rhoR = .125;
    pL = 1; pR = .1;
    
    rho = (rhoL*(x < x0) + rhoR*(x > x0));
    u = 0*x;
    v = 0*x;
    p = (pL*(x < x0) + pR*(x > x0));
end

end