N = 4;
K1D = 8;

CFL = 1; tau = 0;
[err1 entropy rhstest dt Nsteps rho1 rhou1 rhov1 E1] = Euler2D_curved(N,K1D,0,CFL,tau);
dS = abs(entropy-entropy(1));
% dS = entropy;
figure(2)
step = 2*round(.5/CFL);
semilogy(dt*(1:step:Nsteps),dS(1:step:end),'o--');hold on
semilogy(dt*(1:step:Nsteps),abs(rhstest(1:step:end)),'x--');hold on
legend('dS','rhstest')

print_pgf_coordinates(dt*(1:step:Nsteps),dS(1:step:end))
print_pgf_coordinates(dt*(1:step:Nsteps),abs(rhstest(1:step:end)))
return

CFL = .5; tau = 0;
CFL = .25; tau = 1;
[err1 entropy rhstest dt Nsteps rho1 rhou1 rhov1 E1 Vq Pq wJq] = Euler2D_curved(N,K1D,0,CFL,tau);
[err2 entropy rhstest dt Nsteps rho2 rhou2 rhov2 E2] = Euler2D_curved(N,K1D,1,CFL,tau);
[err1,err2]
sqrt(sum(sum(wJq.*((Vq*Pq*(rho1-rho2)).^2  + (Vq*Pq*(rhou1-rhou2)).^2 + (Vq*Pq*(rhov1-rhov2)).^2 + (Vq*Pq*(E1-E2)).^2))))



%%

function [L2err entropy rhstest dt Nsteps rho rhou rhov E Vq2 Pq wJq2] = Euler2D_curved(Nin,K1D,wadgProjEntropyVars,CFL,tauin)

% clear;clear global
Globals2D
global tau

a = 1/4; % warping factor
FinalTime = 2;


if nargin==0
    N = 4;
    K1D = 16;
    wadgProjEntropyVars = 1;
    
    CFL = .5;
    tau = 0;
else
    N = Nin;
    tau = tauin;
end


% Lx = 7.5; Ly = 5; ratiox = 3/4; ratioy = .5;
Lx = 10; Ly = 5; ratiox = 1; ratioy = Ly/Lx;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(round(ratiox*K1D),round(K1D*ratioy));

% ids = abs(abs(VX)-1)>1e-8 & abs(abs(VY)-1)>1e-8;
% VY(ids) = VY(ids) + .05*randn(size(VX(ids)));


VX = VX/max(abs(VX));  VY = VY/max(abs(VY));
VX = (VX+1)*Lx; VY = VY*Ly;

% VX = VX + randn(size(VX));
StartUp2D;
BuildPeriodicMaps2D(max(VX)-min(VX),max(VY)-min(VY));

PlotMesh2D;axis on;return

% plotting nodes
[rp sp] = EquiNodes2D(15); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;
% PlotMesh2D; axis on;return

global M Vq Pq Lq Vfqf Vfq Pfqf VqPq
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ nrJq nsJq wfq
global mapPq
Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
M = Vq'*diag(wq)*Vq;

Pq = M\(Vq'*diag(wq)); % J's cancel out
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

rxJ = rx.*J; sxJ = sx.*J;
ryJ = ry.*J; syJ = sy.*J;
rxJ = Vq*rxJ; sxJ = Vq*sxJ;
ryJ = Vq*ryJ; syJ = Vq*syJ;
J = Vq*J;

[rq1D wq1D] = JacobiGQ(0,0,N);
rfq = [rq1D; -rq1D; -ones(size(rq1D))];
sfq = [-ones(size(rq1D)); rq1D; -rq1D];
wfq = [wq1D; wq1D; wq1D];
Vq1D = Vandermonde1D(N,rq1D)/Vandermonde1D(N,JacobiGL(0,0,N));
% plot(rfq,sfq,'o')
Nfq = length(rq1D);

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
VqPq = Vq*Pq;

global DNr DNs WN
DNr = [Drq .5*VqLq*diag(nrJ);
    -.5*diag(nrJ)*VfPq 0*.5*diag(nrJ)];
DNs = [Dsq .5*VqLq*diag(nsJ);
    -.5*diag(nsJ)*VfPq 0*.5*diag(nsJ)];
WN = diag([wq;wfq]);
% Qr = WN*DNr;
% Qs = WN*DNs;
% keyboard
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
            
            % NOTE - does not work if K1D is too small!!
            if length(p) == 0
%                 keyboard
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
% x0 = 0; y0 = 0; Lx = 1; Ly = 1;
% x = x + Lx*a*cos(1/2*pi*(x-x0)/Lx).*cos(3/2*pi*(y-y0)/Ly);
y = y + Ly*a*sin(pi*(x-x0)/Lx).*cos(1/2*pi*(y-y0)/Ly);
% y(abs(y)<1e-8 & abs(x-10)<1e-8) = y(abs(y)<1e-8 & abs(x-10)<1e-8) + 4*a;
% keyboard


xq = Vq*x; yq = Vq*y;
xp = Vp*x; yp = Vp*y;
xf = Vfq*x;    yf = Vfq*y;

if 0
    rp1D = linspace(-1,1,100)';
    Vp1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,JacobiGL(0,0,N));
    Vfp = kron(eye(Nfaces),Vp1D);
    xfp = Vfp*x(Fmask(:),:);
    yfp = Vfp*y(Fmask(:),:);
    plot(xfp,yfp,'k.')
    hold on
%     plot(x,y,'o')
    axis off
    axis equal
    keyboard
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

% % check compatibility of normals with geofacs
% xf = Vfq*x;
% yf = Vfq*y;

% vv = (nx(mapMq)+nx(mapPq)).^2 + (ny(mapMq)+ny(mapPq)).^2;
% color_line3(xf,yf,vv,vv,'.')
% return

% e = 1;
% [rxJ1 rxJ2] = meshgrid([rxJ(:,e);rxJf(:,e)]);
% [ryJ1 ryJ2] = meshgrid([ryJ(:,e);ryJf(:,e)]);

% Dx = DNr.*(.5*(rxJ1+rxJ2) + .5*(ryJ1+ryJ2));
% diag(wfq.*nrJ.*(rxJf(:,e) + ryJf(:,e)))
% diag(wfq.*nxJ(:,e))
% keyboard

% norm(Dr*Pq*rxJ + Ds*Pq*sxJ,'fro')
% norm(Dr*Pq*ryJ + Ds*Pq*syJ,'fro')

% for i = 1:Nfq*Nfaces*K
%     clf
%     rp1D = linspace(-1,1,100)';
%     Vp1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,JacobiGL(0,0,N));
%     Vfp = kron(eye(Nfaces),Vp1D);
%     xfp = Vfp*x(Fmask(:),:);
%     yfp = Vfp*y(Fmask(:),:);
%     plot(xfp,yfp,'k.')
%     hold on
%
%     plot(xf,yf,'o')
%     plot(xf(i),yf(i),'.','markersize',16)
%     plot(xf(mapPq(i)),yf(mapPq(i)),'x','markersize',16)
%     pause
% end


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


%% test 2 point flux condition

% rhoL = rand+1;
% uL = randn;
% vL = randn;
% EL = 1 + .5*rhoL.*(uL.^2+vL.^2);
% rhoR = 2;
% uR = 3;
% vR = 4;
% ER = 1 + .5*rhoR.*(uR.^2+vR.^2);
%
% V1L = V1(rhoL,rhoL*uL,rhoL*vL,EL);  V1R = V1(rhoR,rhoR*uR,rhoR*vR,ER);
% V2L = V2(rhoL,rhoL*uL,rhoL*vL,EL);  V2R = V2(rhoR,rhoR*uR,rhoR*vR,ER);
% V3L = V3(rhoL,rhoL*uL,rhoL*vL,EL);  V3R = V3(rhoR,rhoR*uR,rhoR*vR,ER);
% V4L = V4(rhoL,rhoL*uL,rhoL*vL,EL);  V4R = V4(rhoR,rhoR*uR,rhoR*vR,ER);
%
%
% val1 = (V1L - V1R)*fxS1(rhoL,uL,vL,EL,rhoR,uR,vR,ER) + ...
%     (V2L - V2R)*fxS2(rhoL,uL,vL,EL,rhoR,uR,vR,ER) + ...
%     (V3L - V3R)*fxS3(rhoL,uL,vL,EL,rhoR,uR,vR,ER) + ...
%     (V4L - V4R)*fxS4(rhoL,uL,vL,EL,rhoR,uR,vR,ER);
%
% val2 = psix(rhoL,rhoL*uL,rhoL*vL,EL) - psix(rhoR,rhoR*uR,rhoR*vR,ER);
%
% norm(val1-val2)
% val1 = (V1L - V1R)*fyS1(rhoL,uL,vL,EL,rhoR,uR,vR,ER) + ...
%     (V2L - V2R)*fyS2(rhoL,uL,vL,EL,rhoR,uR,vR,ER) + ...
%     (V3L - V3R)*fyS3(rhoL,uL,vL,EL,rhoR,uR,vR,ER) + ...
%     (V4L - V4R)*fyS4(rhoL,uL,vL,EL,rhoR,uR,vR,ER);
%
% val2 = psiy(rhoL,rhoL*uL,rhoL*vL,EL) - psiy(rhoR,rhoR*uR,rhoR*vR,ER);
%
% norm(val1-val2)
% keyboard

%% problem params setup

x0 = 0; y0 = 0;

[rhoq uq vq pq] = vortexSolution(xq,yq,0);
rho  = Pq*rhoq;
rhou = Pq*(rhoq.*uq);
rhov = Pq*(rhoq.*vq);
E    = Pq*(pq/(gamma-1) + .5*rhoq.*(uq.^2+vq.^2));

% rho = Pq*(2 + exp(-5^2*(xq).^2));
% rhou = 0*x;
% rhov = 0*x;
% E = 0*x;
% keyboard
% vv = Vp*rho; color_line3(xp,yp,vv,vv,'.'); return

rho = Vq*rho;
rhou = Vq*rhou;
rhov = Vq*rhov;
E = Vq*E;

% test wadg curvi projection property
% e = 1;
% v = Vq*randn(Np,1);
% vJ = VqPq*((VqPq*(v.*J(:,e)))./J(:,e));
% norm(vJ'*diag(wq)*Vq*inv(Vq'*diag(wq./J(:,e))*Vq)*M-(v'*diag(wq.*J(:,e))*Vq))
% keyboard

%%

global wJq
wJq = diag(wq)*(J);

% Runge-Kutta residual storage
res1 = zeros(Nq,K);
res2 = zeros(Nq,K);
res3 = zeros(Nq,K);
res4 = zeros(Nq,K);

% compute time step size
CN = (N+1)*(N+2)/2; % guessing...
CNh = max(CN*max(sJ(:)./Jf(:)));
dt = CFL*1/CNh;

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
        if wadgProjEntropyVars
            q1 = Pq*((VqPq*(V1(rhoq,rhouq,rhovq,Eq).*J))./J);
            q2 = Pq*((VqPq*(V2(rhoq,rhouq,rhovq,Eq).*J))./J);
            q3 = Pq*((VqPq*(V3(rhoq,rhouq,rhovq,Eq).*J))./J);
            q4 = Pq*((VqPq*(V4(rhoq,rhouq,rhovq,Eq).*J))./J);
%             keyboard
        else
            q1 = Pq*V1(rhoq,rhouq,rhovq,Eq);
            q2 = Pq*V2(rhoq,rhouq,rhovq,Eq);
            q3 = Pq*V3(rhoq,rhouq,rhovq,Eq);
            q4 = Pq*V4(rhoq,rhouq,rhovq,Eq);
        end
        
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
        
        % extra LF flux info
        QM{1} = rhoM;       QM{2} = rhouM;
        QM{3} = rhovM;      QM{4} = EM;
        [rhs1 rhs2 rhs3 rhs4]  = RHS2Dsimple(rhoq,uq,vq,Eq,rhoM,uM,vM,EM,QM);
        %         [rhs1 rhs2 rhs3 rhs4]  = RHS2D(rhoq,uq,vq,Eq,rhoM,uM,vM,EM,QM);
        
        if (INTRK==5)
            rhstest(i) = 0;
            for e = 1:K
                %                 rr = wJq(:,e).*(q1q(:,e).*(rhs1(:,e)) + q2q(:,e).*(rhs2(:,e)) + q3q(:,e).*(rhs3(:,e)) + q4q(:,e).*(rhs4(:,e)));
                %                 rhstest(i) = rhstest(i) + sum(rr);               
%                 sum(wJq(:,e).*q1q(:,e).*rhs1(:,e))
                
                r1 = M*((Vq'*diag(wq./J(:,e))*Vq)\(M*Pq*rhs1(:,e)));
                r2 = M*((Vq'*diag(wq./J(:,e))*Vq)\(M*Pq*rhs2(:,e)));
                r3 = M*((Vq'*diag(wq./J(:,e))*Vq)\(M*Pq*rhs3(:,e)));
                r4 = M*((Vq'*diag(wq./J(:,e))*Vq)\(M*Pq*rhs4(:,e)));                
                rhstest(i) = rhstest(i) + sum((q1(:,e)'*(r1) + q2(:,e)'*(r2) + q3(:,e)'*(r3) + q4(:,e)'*(r4)));                               

            end
            

            %             keyboard
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
    
    Sq = -rho.*s(rho,rhou,rhov,E);
    entropy(i) = sum(sum(wJq.*Sq));
    
    if  (mod(i,5)==0 || i==Nsteps)
        clf
        pp = rho;
        vv = real(Vp*Pq*pp);
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f, N = %d, K1D = %d',dt*i,N,K1D))
        %                 view(3)
        drawnow
    end
    
end

[rq2 sq2 wq2] = Cubature2D(Nq+2);
Vq2 = Vandermonde2D(N,rq2,sq2)/V;
xq2 = Vq2*x; yq2 = Vq2*y;
wJq2 = diag(wq2)*(Vq2*Pq*J);
[rhoex uex vex pex] = vortexSolution(xq2,yq2,FinalTime);

rhouex = rhoex.*uex;
rhovex = rhoex.*vex;
% p = (gamma-1)*(E-.5*rho*(u^2+v^2));
Eex = pex/(gamma-1) + .5*rhoex.*(uex.^2+vex.^2);

rhoq = Vq2*Pq*rho;
rhouq = Vq2*Pq*rhou;
rhovq = Vq2*Pq*rhov;
Eq = Vq2*Pq*E;
err = wJq2.*((rhoq-rhoex).^2 + (rhouq-rhouex).^2 + (rhovq-rhovex).^2 + (Eq-Eex).^2);
L2err = sqrt(sum(err(:)));

end


function [rhs1 rhs2 rhs3 rhs4] = RHS2Dsimple(rhoq,uq,vq,Eq,rhoM,uM,vM,EM,QM)

Globals2D;

global M Vq Pq Lq Lqf Vfqf Vfq Pfqf VqPq
% global Drq Dsq Lrq Lsq
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ nrJq nsJq
global mapPq
global fxS1 fyS1 fxS2 fyS2 fxS3 fyS3 fxS4 fyS4
global pfun beta pavg plogmean vnormavg avg
global Drq Dsq VfPq VqLq

global DNr DNs

rhoP = rhoM(mapPq);
uP = uM(mapPq);
vP = vM(mapPq);
EP = EM(mapPq);

betaq = beta(rhoq,uq,vq,Eq);
betafq = beta(rhoM,uM,vM,EM);
betaN = [betaq;betafq];

% Lax-Friedrichs flux
global gamma tau
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

f1 = nxJ.*fxS1(rhoM,uM,vM,EM,rhoM,uM,vM,EM) + nyJ.*fyS1(rhoM,uM,vM,EM,rhoM,uM,vM,EM);
f2 = nxJ.*fxS2(rhoM,uM,vM,EM,rhoM,uM,vM,EM) + nyJ.*fyS2(rhoM,uM,vM,EM,rhoM,uM,vM,EM);
f3 = nxJ.*fxS3(rhoM,uM,vM,EM,rhoM,uM,vM,EM) + nyJ.*fyS3(rhoM,uM,vM,EM,rhoM,uM,vM,EM);
f4 = nxJ.*fxS4(rhoM,uM,vM,EM,rhoM,uM,vM,EM) + nyJ.*fyS4(rhoM,uM,vM,EM,rhoM,uM,vM,EM);

rho = [rhoq; rhoM];
u = [uq; uM];
v = [vq; vM];
E = [Eq; EM];
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
    
    [rxJ1, rxJ2] = meshgrid([rxJ(:,e);rxJf(:,e)]);
    [sxJ1, sxJ2] = meshgrid([sxJ(:,e);sxJf(:,e)]);
    [ryJ1, ryJ2] = meshgrid([ryJ(:,e);ryJf(:,e)]);
    [syJ1, syJ2] = meshgrid([syJ(:,e);syJf(:,e)]);
    rxJK = avg(rxJ1,rxJ2);  sxJK = avg(sxJ1,sxJ2);
    ryJK = avg(ryJ1,ryJ2);  syJK = avg(syJ1,syJ2);
    
    Dx = DNr.*rxJK + DNs.*sxJK;
    Dy = DNr.*ryJK + DNs.*syJK;
    
    divF1(:,e) = sum(Dx.*FxS1,2) + sum(Dy.*FyS1,2);
    divF2(:,e) = sum(Dx.*FxS2,2) + sum(Dy.*FyS2,2);
    divF3(:,e) = sum(Dx.*FxS3,2) + sum(Dy.*FyS3,2);
    divF4(:,e) = sum(Dx.*FxS4,2) + sum(Dy.*FyS4,2);

end

rhs1 = 2*[VqPq VqLq]*divF1 + VqLq*(fSf1-0*f1);
rhs2 = 2*[VqPq VqLq]*divF2 + VqLq*(fSf2-0*f2);
rhs3 = 2*[VqPq VqLq]*divF3 + VqLq*(fSf3-0*f3);
rhs4 = 2*[VqPq VqLq]*divF4 + VqLq*(fSf4-0*f4);

% apply wadg
rhs1 = -VqPq*(rhs1./J);
rhs2 = -VqPq*(rhs2./J);
rhs3 = -VqPq*(rhs3./J);
rhs4 = -VqPq*(rhs4./J);


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

if 1
    % pulse condition
    x0 = 5;
    rho = 2 + (abs(x-x0) < 5);
    u = 0*rho;
    v = 0*rho;
    p = rho.^gamma;
end

end