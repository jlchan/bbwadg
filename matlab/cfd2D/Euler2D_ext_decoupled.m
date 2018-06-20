clear
clear -globals
Globals2D

a = 0/8; % warping factor
% FinalTime = 10;

N = 4;
K1D = 8;
FinalTime = 1;

wadgProjEntropyVars = abs(a)>1e-8;
CFL = .25;
global tau
tau = 0;

% Lx = 7.5; Ly = 5; ratiox = 3/4; ratioy = .5;
Lx = 10; Ly = 5; ratiox = 1; ratioy = Ly/Lx;
%[Nv, VX, VY, K, EToV] = unif_tri_mesh(round(ratiox*K1D),round(K1D*ratioy));
[Nv, VX, VY, K, EToV] = QuadMesh2D(round(ratiox*K1D),round(K1D*ratioy));

% ids = abs(abs(VX)-1)>1e-8 & abs(abs(VY)-1)>1e-8;
% VY(ids) = VY(ids) + .05*randn(size(VX(ids)));


VX = VX/max(abs(VX));  VY = VY/max(abs(VY));
VX = (VX+1)*Lx; VY = VY*Ly;

% VX = VX + randn(size(VX));
StartUp2D;
BuildPeriodicMaps2D(max(VX)-min(VX),max(VY)-min(VY));

% PlotMesh2D;axis on;return

% plotting nodes
[rp sp] = EquiNodes2D(10); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;
% PlotMesh2D; axis on;return

global M Vq Pq Lq Vfqf Vfq Pfqf VqPq VqLq
global nrJ nsJ nrJq nsJq
% Nq = 2*N+1;
% [rq sq wq] = Cubature2D(Nq); % integrate u*v*c
% [rq sq wq] = QNodes2D(N); if length(rq)~=Np;  keyboard; end
[rq1D wq1D] = JacobiGQ(0,0,N);
[rq sq] = meshgrid(rq1D);
rq = rq(:); sq = sq(:);
[wr ws] = meshgrid(wq1D);
wq = wr(:).*ws(:);
Nq = length(wq);

Vq = Vandermonde2D(N,rq,sq)/V;
M = Vq'*diag(wq)*Vq;

Pq = M\(Vq'*diag(wq)); % J's cancel out
xq = Vq*x; yq = Vq*y;

% face nodes
% [rq1D wq1D] = JacobiGQ(0,0,N+1);
rq1D = [-1 + (1+rq1D)/2; (1+rq1D)/2];
wq1D = .5*[wq1D;wq1D];
Nfq = length(rq1D);
e = ones(size(rq1D));
rfq = [rq1D; e; -rq1D; -e];
sfq = [-e; rq1D; e; -rq1D];
wfq = [wq1D; wq1D; wq1D; wq1D];
V1D = Vandermonde1D(N,JacobiGL(0,0,N));
Vq1D = Vandermonde1D(N,rq1D)/V1D;
Vfqf = kron(eye(Nfaces),Vq1D);
Vfq = Vandermonde2D(N,rfq,sfq)/V;
Mf = Vfq'*diag(wfq)*Vfq;
Lq = M\(Vfq'*diag(wfq));

% multi-layer decoupling - a "middle" layer
[rq1Dmid wq1Dmid] = JacobiGQ(0,0,N);
Nfqmid = length(rq1Dmid);
emid = ones(size(rq1Dmid));
rfqmid = [rq1Dmid; emid; -rq1Dmid; -emid];
sfqmid = [-emid; rq1Dmid; emid; -rq1Dmid];
wfqmid = [wq1Dmid; wq1Dmid; wq1Dmid; wq1Dmid];
Vq1Dmid = Vandermonde1D(N,rq1Dmid)/V1D;
VfPfmid = Vq1D*((Vq1Dmid'*diag(wq1Dmid)*Vq1Dmid)\(Vq1Dmid'*diag(wq1Dmid)));
VfPfmid = kron(eye(Nfaces),VfPfmid);
VfmidPf = Vq1Dmid*((Vq1D'*diag(wq1D)*Vq1D)\(Vq1D'*diag(wq1D)));
VfmidPf = kron(eye(Nfaces),VfmidPf);

global Cf Lqmid VqLqmid
Cf = [zeros(length(rfqmid)) VfmidPf; -VfPfmid zeros(length(rfq))];
Sf = diag([wfqmid;wfq])*Cf; % skew version 
Vfqmid = Vandermonde2D(N,rfqmid,sfqmid)/V;
Vfqfmid = kron(eye(Nfaces),Vq1Dmid);
Lqmid = M\(Vfqmid'*diag(wfqmid));
norm(Sf+Sf','fro')

nrJ = [0*e; e; 0*e; -e]; % sJ = 2 for all faces,
nsJ = [-e; 0*e; e; 0*e];
nrJmid = [0*emid; emid; 0*emid; -emid]; % sJ = 2 for all faces,
nsJmid = [-emid; 0*emid; emid; 0*emid];
nrJq = repmat(nrJmid',Nq,1);
nsJq = repmat(nsJmid',Nq,1);

% flux differencing operators
VqPq = Vq*Pq;
VfmidPq = (Vfqmid*Pq);
VqLqmid = Vq*Lqmid;
VqLq = Vq*Lq;

global DNr DNs WN
wN = [wq; wfqmid];
WN = diag(wN);

Qr = diag(wq)*Vq*Dr*Pq;
Qs = diag(wq)*Vq*Ds*Pq;

QNrskew = .5*[Qr-Qr' diag(wq)*VqLqmid*diag(nrJmid);
    -diag(wfqmid.*nrJmid)*VfmidPq zeros(length(nsJmid))];
QNsskew = .5*[Qs-Qs' diag(wq)*VqLqmid*diag(nsJmid);
    -diag(wfqmid.*nsJmid)*VfmidPq zeros(length(nsJmid))];

% if 0
%     wN = [wq; wfqmid; wfq];
%     Z13 = zeros(size(Qr,1),size(VfmidPf,2));
%     Z22 = zeros(length(nrJmid));
%     Z33 = zeros(length(nrJ));
%     QNrskew = .5*[Qr-Qr' diag(wq)*VqLq*diag(nrJmid) Z13;
%         -diag(wfqmid.*nrJmid)*VfPq Z22 VfmidPf*diag(wfq.*nrJ);
%         Z13' -diag(wfq)*VfPfmid*diag(nrJmid) Z33 ];
%     
%     QNsskew = .5*[Qs-Qs' diag(wq)*VqLq*diag(nsJmid) Z13;
%         -diag(wfqmid.*nsJmid)*VfPq Z22 diag(wfqmid.*nsJmid)*VfmidPf;
%         Z13' -diag(wfq)*VfPfmid*diag(nsJmid) Z33];
%     
%     % define face interpolation points as xf, yf. ASSUMES DISTINCT POINTS
%     Vfq = [Vfqmid;Vfq];
%     VqLq = Vq*(M\(Vfq'*diag([wfqmid;wfq])));
%     nrJ = [nrJmid;nrJ];
%     nsJ = [nsJmid;nsJ];
%     midmask = [0*wfqmid;ones(size(wfq))];
% end

% make skew symmetric diff matrices
DNr = diag(1./wN)*QNrskew;
DNs = diag(1./wN)*QNsskew;
% return


%% make quadrature face maps

global mapPq

xf = Vfq*x;
yf = Vfq*y;

% plot(xf,yf,'o');return

Nfq = size(Vfq,1)/Nfaces; % redefine Nfq for enriched space
mapMq = reshape(1:length(xf(:)),Nfq*Nfaces,K);
mapPq = mapMq;
% keyboard

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


%% make curvilinear mesh

x0 = Lx; y0 = 0;
% x0 = 0; y0 = 0; Lx = 1; Ly = 1;
x = x + Lx*a*cos(1/2*pi*(x-x0)/Lx).*cos(3/2*pi*(y-y0)/Ly);
y = y + Ly*a*sin(pi*(x-x0)/Lx).*cos(1/2*pi*(y-y0)/Ly);

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
    plot(x,y,'o')
    L2err = nan;
    return
end

global rxJ sxJ ryJ syJ 
global rxJf sxJf ryJf syJf
global nxJ nyJ 

rxJ = zeros(Nq,K); sxJ = zeros(Nq,K);
ryJ = zeros(Nq,K); syJ = zeros(Nq,K);
J = zeros(Nq,K);

NfqNfaces = size(Vfq,1);
rxJf = zeros(NfqNfaces,K); sxJf = zeros(NfqNfaces,K);
ryJf = zeros(NfqNfaces,K); syJf = zeros(NfqNfaces,K);
Jf = zeros(NfqNfaces,K);

NfqNfaces = size(Vfqmid,1);
rxJfmid = zeros(NfqNfaces,K); sxJfmid = zeros(NfqNfaces,K);
ryJfmid = zeros(NfqNfaces,K); syJfmid = zeros(NfqNfaces,K);
Jfmid = zeros(NfqNfaces,K);

for e = 1:K
    [rxk,sxk,ryk,syk,Jk] = GeometricFactors2D(x(:,e),y(:,e),Vq*Dr,Vq*Ds);
    rxJ(:,e) = rxk.*Jk;    sxJ(:,e) = sxk.*Jk;
    ryJ(:,e) = ryk.*Jk;    syJ(:,e) = syk.*Jk;
    J(:,e) = Jk;       
    
    [rxk,sxk,ryk,syk,Jk] = GeometricFactors2D(x(:,e),y(:,e),Vfq*Dr,Vfq*Ds);
    rxJf(:,e) = rxk.*Jk;   sxJf(:,e) = sxk.*Jk;
    ryJf(:,e) = ryk.*Jk;   syJf(:,e) = syk.*Jk;
    Jf(:,e) = Jk;
    
    [rxk,sxk,ryk,syk,Jk] = GeometricFactors2D(x(:,e),y(:,e),Vfqmid*Dr,Vfqmid*Ds);
    rxJfmid(:,e) = rxk.*Jk;   sxJfmid(:,e) = sxk.*Jk;
    ryJfmid(:,e) = ryk.*Jk;   syJfmid(:,e) = syk.*Jk;
    Jfmid(:,e) = Jk;
end

global nxJmid nyJmid
nxJmid = rxJfmid.*nrJmid + sxJfmid.*nsJmid;
nyJmid = ryJfmid.*nrJmid + syJfmid.*nsJmid;

nxJ = rxJf.*nrJ + sxJf.*nsJ;
nyJ = ryJf.*nrJ + syJf.*nsJ;

rxJf = rxJfmid; sxJf = sxJfmid;
ryJf = ryJfmid; syJf = syJfmid;

nx = nxJ./Jf;
ny = nyJ./Jf;
sJ = sqrt(nx.^2 + ny.^2);
nx = nx./sJ; ny = ny./sJ;
sJ = sJ.*Jf;

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
        
        % evaluate at quad/surface points
        q1q = Vq*q1;  
        q2q = Vq*q2;  
        q3q = Vq*q3;  
        q4q = Vq*q4;  
        
        rhoq  = U1(q1q,q2q,q3q,q4q); 
        rhouq = U2(q1q,q2q,q3q,q4q); 
        rhovq = U3(q1q,q2q,q3q,q4q); 
        Eq    = U4(q1q,q2q,q3q,q4q); 
        uq = rhouq./rhoq; 
        vq = rhovq./rhoq; 
        
        % intermediate "middle layer" points
        q1M = Vfqmid*q1;
        q2M = Vfqmid*q2;
        q3M = Vfqmid*q3;
        q4M = Vfqmid*q4;                
        rhof  = U1(q1M,q2M,q3M,q4M);
        rhouf = U2(q1M,q2M,q3M,q4M);
        rhovf = U3(q1M,q2M,q3M,q4M);
        Ef    = U4(q1M,q2M,q3M,q4M);        
        uf = rhouf./rhof;
        vf = rhovf./rhof;
        
        % actual points
        q1M = Vfq*q1;
        q2M = Vfq*q2;
        q3M = Vfq*q3;
        q4M = Vfq*q4;                
        rhoM  = U1(q1M,q2M,q3M,q4M);
        rhouM = U2(q1M,q2M,q3M,q4M);
        rhovM = U3(q1M,q2M,q3M,q4M);
        EM    = U4(q1M,q2M,q3M,q4M);        
        uM = rhouM./rhoM;
        vM = rhovM./rhoM;
        
        % extra LF flux info
        [rhs1 rhs2 rhs3 rhs4]  = RHS2Dsimple(rhoq,uq,vq,Eq,rhof,uf,vf,Ef,rhoM,uM,vM,EM);
        
        if (INTRK==5)
            rhstest(i) = 0;
            for e = 1:K
                %                 r1 = M*((Vq'*diag(wq./J(:,e))*Vq)\(M*Pq*rhs1(:,e)));
                %                 r2 = M*((Vq'*diag(wq./J(:,e))*Vq)\(M*Pq*rhs2(:,e)));
                %                 r3 = M*((Vq'*diag(wq./J(:,e))*Vq)\(M*Pq*rhs3(:,e)));
                %                 r4 = M*((Vq'*diag(wq./J(:,e))*Vq)\(M*Pq*rhs4(:,e)));
                %                 rhstest(i) = rhstest(i) + sum((q1(:,e)'*(r1) + q2(:,e)'*(r2) + q3(:,e)'*(r3) + q4(:,e)'*(r4)));
                q1e = J(:,e).*wq.*(Vq*q1(:,e));
                q2e = J(:,e).*wq.*(Vq*q2(:,e));
                q3e = J(:,e).*wq.*(Vq*q3(:,e));
                q4e = J(:,e).*wq.*(Vq*q4(:,e));
                r1 = rhs1(:,e);
                r2 = rhs2(:,e);
                r3 = rhs3(:,e);
                r4 = rhs4(:,e);
                rhstest(i) = rhstest(i) + (q1e'*r1 + q2e'*r2 + q3e'*r3 + q4e'*r4);
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
L2err = sqrt(sum(err(:)))

dS = abs(entropy-entropy(1));
dS = entropy + max(abs(entropy));
figure(2)
semilogy(dt*(1:Nsteps),dS,'o--');hold on
semilogy(dt*(1:Nsteps),abs(rhstest),'x--');hold on
legend('dS','rhstest')


function [rhs1 rhs2 rhs3 rhs4] = RHS2Dsimple(rhoq,uq,vq,Eq,rhof,uf,vf,Ef,rhoM,uM,vM,EM)

Globals2D;

global M Vq Pq Lq Lqf Vfqf Vfq VfPq VqPq VqLq
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ nrJq nsJq
global mapPq
global fxS1 fyS1 fxS2 fyS2 fxS3 fyS3 fxS4 fyS4
global pfun beta pavg plogmean vnormavg avg

global DNr DNs

rhoP = rhoM(mapPq);
uP = uM(mapPq);
vP = vM(mapPq);
EP = EM(mapPq);

betaq = beta(rhoq,uq,vq,Eq);
betafq = beta(rhof,uf,vf,Ef);
betaM = beta(rhoM,uM,vM,EM);

% Lax-Friedrichs flux
global gamma tau
unorm2 = (uM.^2+vM.^2);
pM = (gamma-1)*(EM - .5*rhoM.*unorm2);
cvel = sqrt(gamma*pM./rhoM);
lam = sqrt(unorm2)+cvel;
LFc = max(lam(mapPq),lam);

% rhoUn = (rhoP.*uP-rhoM.*uM).*nx + (rhoP.*vP-rhoM.*vM).*ny;
dQ1 = rhoP-rhoM; % QP{1}-QM{1};
dQ2 = rhoP.*uP - rhoM.*uM; % QP{2}-QM{2};
dQ3 = rhoP.*vP - rhoM.*vM; % QP{3}-QM{3};
dQ4 = EP - EM; % QP{4}-QM{4};
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

rho = [rhoq; rhof];
u = [uq; uf];
v = [vq; vf];
betaN = [betaq;betafq];

divF1 = zeros(size(DNr,1),K);
divF2 = zeros(size(DNr,1),K);
divF3 = zeros(size(DNr,1),K);
divF4 = zeros(size(DNr,1),K);
for e = 1:K
    
    [rhox, rhoy] = meshgrid(rho(:,e));
    [ux, uy] = meshgrid(u(:,e));
    [vx, vy] = meshgrid(v(:,e));
    [betax, betay] = meshgrid(betaN(:,e));
    
    % optimized evaluations
    [FxS1 FxS2 FxS3 FxS4 FyS1 FyS2 FyS3 FyS4] = euler_flux(rhox,rhoy,ux,uy,vx,vy,betax,betay);
        
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

global Cf VqLqmid nxJmid nyJmid

VqPqLq = 2*[VqPq VqLqmid];
rhs1 = VqPqLq*divF1 + VqLq*(fSf1);
rhs2 = VqPqLq*divF2 + VqLq*(fSf2);
rhs3 = VqPqLq*divF3 + VqLq*(fSf3);
rhs4 = VqPqLq*divF4 + VqLq*(fSf4);

for e = 1:K
    [rhox rhoy] = meshgrid([rhof(:,e);rhoM(:,e)]);
    [ux uy] = meshgrid([uf(:,e);uM(:,e)]);
    [vx vy] = meshgrid([vf(:,e);vM(:,e)]);
    [betax betay] = meshgrid([betafq(:,e);betaM(:,e)]);       
    
    [nx1 nx2] = meshgrid([nxJmid(:,e);nxJ(:,e)]);
    [ny1 ny2] = meshgrid([nyJmid(:,e);nyJ(:,e)]);
    nx = avg(nx1,nx2); ny = avg(ny1,ny2);
    [FxS1 FxS2 FxS3 FxS4 FyS1 FyS2 FyS3 FyS4] = euler_flux(rhox,rhoy,ux,uy,vx,vy,betax,betay);
    
    FnS1 = FxS1.*nx + FyS1.*ny;
    FnS2 = FxS2.*nx + FyS2.*ny;
    FnS3 = FxS3.*nx + FyS3.*ny;
    FnS4 = FxS4.*nx + FyS4.*ny;
    VqLqLq = [VqLqmid VqLq];

    rhs1(:,e) = rhs1(:,e) + VqLqLq*sum(Cf.*FnS1,2);
    rhs2(:,e) = rhs2(:,e) + VqLqLq*sum(Cf.*FnS2,2);
    rhs3(:,e) = rhs3(:,e) + VqLqLq*sum(Cf.*FnS3,2);
    rhs4(:,e) = rhs4(:,e) + VqLqLq*sum(Cf.*FnS4,2);
end

% apply wadg at the end
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

function [FxS1 FxS2 FxS3 FxS4 FyS1 FyS2 FyS3 FyS4] = euler_flux(rhox,rhoy,ux,uy,vx,vy,betax,betay)

global avg gamma
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

end