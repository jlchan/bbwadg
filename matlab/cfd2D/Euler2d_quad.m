clear
clear -globals
Globals2D

a = .125; % warping factor

N = 2;
K1D = 4;
FinalTime = 1.0;

wadgProjEntropyVars = abs(a)>1e-8;
CFL = .5;
global tau
tau = 0;
plotMesh = 1;

% Lx = 7.5; Ly = 5; ratiox = 3/4; ratioy = .5;
Lx = 10; Ly = 5; ratiox = 1; ratioy = Ly/Lx;
%[Nv, VX, VY, K, EToV] = unif_tri_mesh(round(ratiox*K1D),round(K1D*ratioy));
[Nv, VX, VY, K, EToV] = QuadMesh2D(2*K1D,K1D);

% ids = abs(abs(VX)-1)>1e-8 & abs(abs(VY)-1)>1e-8;
% VY(ids) = VY(ids) + .05*randn(size(VX(ids)));

VX = VX/max(abs(VX));  VY = VY/max(abs(VY));
VX = (VX+1)*Lx; VY = VY*Ly;

% VX = VX + .75*randn(size(VX));
% VY = VY + .75*randn(size(VX));
StartUp2D;
BuildPeriodicMaps2D(max(VX)-min(VX),max(VY)-min(VY));


% plotting nodes
[rp sp] = EquiNodes2D(15); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;
% PlotMesh2D; axis on;return

global M Vq Pq Lq Vfqf Vfq Pfqf VqPq VqLq
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ wfq
global mapPq

% vol nodes
[rq1D wq1D] = JacobiGL(0,0,N); 
[rq sq] = meshgrid(rq1D);
rq = rq(:); sq = sq(:);
[wr ws] = meshgrid(wq1D); 
wq = wr(:).*ws(:);

Nq = length(rq);

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

% face nodes
[rq1D wq1D] = JacobiGQ(0,0,N); % 2N-1
% rq1D = rq1D*(1-1e-10);
% rq1D = [-1+(1+rq1D)/2; (1+rq1D)/2];
% wq1D = [wq1D; wq1D]/2;
% rq1D = rq1D*(1-1e-9);
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

nx = Vfqf*nx;
ny = Vfqf*ny;
sJ = Vfqf*sJ;
nxJ = (nx.*sJ);
nyJ = (ny.*sJ);
Fscale = Vfqf*Fscale;

nrJ = [0*e; e; 0*e; -e]; % sJ = 2 for all faces, 
nsJ = [-e; 0*e; e; 0*e];
% quiver(rfq,sfq,nrJ,nsJ);return

% flux differencing operators

Drq = (Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vfq*Pq);
Dsq = (Vq*Ds*Pq - .5*Vq*Lq*diag(nsJ)*Vfq*Pq);
VfPq = (Vfq*Pq);
VqLq = Vq*Lq;
VqPq = Vq*Pq;

global DNr DNs WN 
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

% keyboard
% return

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

% Qr = diag(wq)*(Vq*Dr*Pq);
% norm(Qr+Qr' - Vfq'*diag(nrJ.*wfq)*Vfq)


%% make quadrature face maps

xf = Vfq*x;
yf = Vfq*y;

% mapMq = reshape(mapM,Nfp*Nfaces,K);
% mapPq = reshape(mapP,Nfp*Nfaces,K);
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


%% make curvilinear mesh 

x0 = Lx; y0 = 0;
% x0 = 0; y0 = 0; Lx = 1; Ly = 1;
dx = cos(1/2*pi*(x-x0)/Lx).*cos(3/2*pi*(y-y0)/Ly);
x = x + Lx*a*dx;
dy = sin(4/2*pi*(x-x0)/Lx).*cos(1/2*pi*(y-y0)/Ly);
y = y + Ly*a*dy;

xq = Vq*x; yq = Vq*y;
xp = Vp*x; yp = Vp*y;
xf = Vfq*x;    yf = Vfq*y;

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
% e = 1;
% DNx = .5*(diag([rxJ(:,e); rxJf(:,e)])*DNr + DNr*diag([rxJ(:,e); rxJf(:,e)]) ...
%      + diag([sxJ(:,e); sxJf(:,e)])*DNs + DNs*diag([sxJ(:,e); sxJf(:,e)]));
% BNr = [zeros(Nq) zeros(Nq,Nfq*Nfaces);
%     zeros(Nfq*Nfaces,Nq) diag(wfq.*nrJ)];
% BNs = [zeros(Nq) zeros(Nq,Nfq*Nfaces);
%     zeros(Nfq*Nfaces,Nq) diag(wfq.*nsJ)];
% 
% BNx = diag([rxJ(:,e); rxJf(:,e)])*BNr + diag([sxJ(:,e); sxJf(:,e)])*BNs;
% QNx = WN*DNx;
% 
% u = randn(size(DNx,2),1);
% e = ones(size(DNx,2),1);
% e'*QNx*u
% e'*(BNx-QNx')*u

% keyboard

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
CN = (N+1)*(N+2)/2; % for quads
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
        
        %         % project to entropy variables for curvilinear
        %         if wadgProjEntropyVars
        %             q1 = Pq*((VqPq*(V1(rhoq,rhouq,rhovq,Eq).*J))./J);
        %             q2 = Pq*((VqPq*(V2(rhoq,rhouq,rhovq,Eq).*J))./J);
        %             q3 = Pq*((VqPq*(V3(rhoq,rhouq,rhovq,Eq).*J))./J);
        %             q4 = Pq*((VqPq*(V4(rhoq,rhouq,rhovq,Eq).*J))./J);
        %         else
        %             q1 = Pq*V1(rhoq,rhouq,rhovq,Eq);
        %             q2 = Pq*V2(rhoq,rhouq,rhovq,Eq);
        %             q3 = Pq*V3(rhoq,rhouq,rhovq,Eq);
        %             q4 = Pq*V4(rhoq,rhouq,rhovq,Eq);
        %         end
        
        % evaluate at quad/surface points
        %         q1q = Vq*q1;
        %         q2q = Vq*q2;
        %         q3q = Vq*q3;
        %         q4q = Vq*q4;
        q1q = V1(rhoq,rhouq,rhovq,Eq);
        q2q = V2(rhoq,rhouq,rhovq,Eq);
        q3q = V3(rhoq,rhouq,rhovq,Eq);
        q4q = V4(rhoq,rhouq,rhovq,Eq);
        
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
        QM{1} = rhoM;       QM{2} = rhouM;
        QM{3} = rhovM;      QM{4} = EM;
        [rhs1 rhs2 rhs3 rhs4]  = RHS2Dsimple(rhoq,uq,vq,Eq,rhoM,uM,vM,EM,QM);        
        
        if (INTRK==5)
            rhstest(i) = sum(sum(wJq.*(q1q.*rhs1 + q2q.*rhs2 + q3q.*rhs3 + q4q.*rhs4)));
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

[rq2 sq2 wq2] = Cubature2D(Nq+1);
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

% return
dS = abs(entropy-entropy(1));
% dS = entropy + max(abs(entropy));
figure(2)
semilogy(dt*(1:Nsteps),dS,'o--');hold on
semilogy(dt*(1:Nsteps),abs(rhstest),'x--');hold on
legend('dS','rhstest')


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

betaq = beta(rhoq,uq,vq,Eq);
betafq = beta(rhoM,uM,vM,EM);

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

% f1 = nxJ.*fxS1(rhoM,uM,vM,EM,rhoM,uM,vM,EM) + nyJ.*fyS1(rhoM,uM,vM,EM,rhoM,uM,vM,EM);
% f2 = nxJ.*fxS2(rhoM,uM,vM,EM,rhoM,uM,vM,EM) + nyJ.*fyS2(rhoM,uM,vM,EM,rhoM,uM,vM,EM);
% f3 = nxJ.*fxS3(rhoM,uM,vM,EM,rhoM,uM,vM,EM) + nyJ.*fyS3(rhoM,uM,vM,EM,rhoM,uM,vM,EM);
% f4 = nxJ.*fxS4(rhoM,uM,vM,EM,rhoM,uM,vM,EM) + nyJ.*fyS4(rhoM,uM,vM,EM,rhoM,uM,vM,EM);

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
    
    % apply correction term for faces
    
%     divF1(:,e) = divF1(:,e) + sum(Sfx.*Fx
end

PN = 2*[VqPq VqLq];
rhs1 = PN*divF1 + VqLq*(fSf1);
rhs2 = PN*divF2 + VqLq*(fSf2);
rhs3 = PN*divF3 + VqLq*(fSf3);
rhs4 = PN*divF4 + VqLq*(fSf4);

% apply wadg at the end
rhs1 = -(rhs1./J);
rhs2 = -(rhs2./J);
rhs3 = -(rhs3./J);
rhs4 = -(rhs4./J);

% rhs1 = -VqPq*(rhs1./J);
% rhs2 = -VqPq*(rhs2./J);
% rhs3 = -VqPq*(rhs3./J);
% rhs4 = -VqPq*(rhs4./J);

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