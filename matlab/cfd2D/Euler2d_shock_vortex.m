
Globals2D
global tau

a = 0/4; % warping factor
FinalTime = .7;
% FinalTime = .3;

N = 3;
K1D = 48; 
wadgProjEntropyVars = abs(a) > 1e-10;

CFL = .5;
tau = 1;

% Lx = 7.5; Ly = 5; ratiox = 3/4; ratioy = .5;
Lx = 1.5; Ly = 1; ratiox = 1; ratioy = Ly/Lx;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(round(ratiox*K1D),round(K1D*ratioy));
VX = VX/max(abs(VX));  VY = VY/max(abs(VY));
VX = (VX+1)/2*Lx; VY = (VY+1)/2*Ly;

StartUp2D;
BuildPeriodicMaps2D(max(VX)-min(VX), max(VY)-min(VY));

% plotting nodes
Nplot = 15;
[rp sp] = EquiNodes2D(Nplot); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;
% PlotMesh2D; axis on;return

global xq yq 
global M Vq Pq Lq Vfqf Vfq Pfqf VqPq
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ nrJq nsJq wfq
global mapPq
global wq
Nq = 2*N;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
M = Vq'*diag(wq)*Vq;

Pq = M\(Vq'*diag(wq)); % J's cancel out
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

[rq1D wq1D] = JacobiGQ(0,0,N);
rfq = [rq1D; -rq1D; -ones(size(rq1D))];
sfq = [-ones(size(rq1D)); rq1D; -rq1D];
wfq = [wq1D; sqrt(8)*wq1D; wq1D];
Vq1D = Vandermonde1D(N,rq1D)/Vandermonde1D(N,JacobiGL(0,0,N));
% plot(rfq,sfq,'o')
Nfq = length(rq1D);

Vfq = Vandermonde2D(N,rfq,sfq)/V;
Vfqf = kron(eye(3),Vq1D);
Mf = Vfq'*diag(wfq)*Vfq;
Lq = M\(Vfq'*diag(wfq));

Pq1D = (Vq1D'*diag(wq1D)*Vq1D) \ (Vq1D'*diag(wq1D));
Pfqf = kron(eye(3),Pq1D);

nrJ = [-zeros(size(rq1D)); ones(size(rq1D)); -ones(size(rq1D)); ]; 
nsJ = [-ones(size(rq1D)); ones(size(rq1D)); -zeros(size(rq1D)); ];

Nq = length(rq);
nrJq = repmat(nrJ',Nq,1);
nsJq = repmat(nsJ',Nq,1);

% recompute geofacs
rxJ = rx.*J; sxJ = sx.*J;
ryJ = ry.*J; syJ = sy.*J;
rxJ = Vq*rxJ; sxJ = Vq*sxJ;
ryJ = Vq*ryJ; syJ = Vq*syJ;
J = Vq*J;
nx = Vfqf*nx;
ny = Vfqf*ny;
sJ = Vfqf*sJ;
nxJ = (nx.*sJ);
nyJ = (ny.*sJ);
Fscale = Vfqf*Fscale;

% flux differencing operators
global Drq Dsq VfPq VqLq
Drq = (Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vfq*Pq);
Dsq = (Vq*Ds*Pq - .5*Vq*Lq*diag(nsJ)*Vfq*Pq);
VfPq = (Vfq*Pq);
VqLq = Vq*Lq;
VqPq = Vq*Pq;

global DNr DNs WN

% skew formulation matrices
DNr = [Drq .5*VqLq*diag(nrJ);
    -.5*diag(nrJ)*VfPq zeros(length(nrJ))];
DNs = [Dsq .5*VqLq*diag(nsJ);
    -.5*diag(nsJ)*VfPq zeros(length(nrJ))];
WN = diag([wq;wfq]);

QNr = WN*DNr;
QNs = WN*DNs;

% keyboard
% PlotMesh2D;axis on; grid on;return
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

global mapBwall
mapBwall = find(abs(yf-max(y(:)))<1e-8 | abs(yf-min(y(:)))<1e-8);


%% filter

% F = Filter2D(N,.5,16);

%% make curvilinear mesh (still unstable?)

x0 = Lx; y0 = 0;
% x0 = 0; y0 = 0; Lx = 1; Ly = 1;
x = x + .5*Lx*a*cos(1/2*pi*(x-x0)/Lx).*cos(3/2*pi*(y-y0)/Ly);
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
%     keyboard
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

% to ensure weight-adjusted DG is conservative
J = VqPq*J;


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
% keyboard
rho  = Pq*rhoq;
rhou = Pq*(rhoq.*uq);
rhov = Pq*(rhoq.*vq);
E    = Pq*(pq/(gamma-1) + .5*rhoq.*(uq.^2+vq.^2));

% rho = 2+0*Pq*(2 + exp(-5^2*(xq).^2));
% rhou = 0+0*x;
% rhov = 0+0*x;
% pq = (2 + 0*exp(-5^2*(xq).^2)).^gamma;
% E    = Pq*(pq/(gamma-1));
% keyboard
% vv = Vp*rho; color_line3(xp,yp,vv,vv,'.'); return

rho = Vq*rho;
rhou = Vq*rhou;
rhov = Vq*rhov;
E = Vq*E;

%% boundary conditions

global inflow wall 
inflow = find(abs(xf)<1e-8);
wall = find(abs(yf)<1e-8 | abs(yf-1)<1e-8);

% PlotMesh2D;axis on;plot(xf(wall),yf(wall),'o');hold on;return
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

clear Uavg UavgWADG rhsavg rhsavgWADG
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
        [rhs1 rhs2 rhs3 rhs4]  = RHS2Dsimple(rhoq,uq,vq,Eq,rhoM,uM,vM,EM,QM,i*dt);
        %         [rhs1 rhs2 rhs3 rhs4] = RHS2D(rhoq,uq,vq,Eq,rhoM,uM,vM,EM,QM);                
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        res4 = rk4a(INTRK)*res4 + dt*rhs4;
        
        rho  = rho  + rk4b(INTRK)*res1;
        rhou = rhou + rk4b(INTRK)*res2;
        rhov = rhov + rk4b(INTRK)*res3;
        E    = E    + rk4b(INTRK)*res4;
        
    end;    
    
    %     rhstest(i) = 0;
    %     rhsavgWADG(i) = 0;
    %     UavgWADG(i) = 0;
    %     for e = 1:K
    %         r1 = M*((Vq'*diag(wq./J(:,e))*Vq)\(M*Pq*rhs1(:,e)));
    %         r2 = M*((Vq'*diag(wq./J(:,e))*Vq)\(M*Pq*rhs2(:,e)));
    %         r3 = M*((Vq'*diag(wq./J(:,e))*Vq)\(M*Pq*rhs3(:,e)));
    %         r4 = M*((Vq'*diag(wq./J(:,e))*Vq)\(M*Pq*rhs4(:,e)));
    %         rhstest(i) = rhstest(i) + sum((q1(:,e)'*(r1) + q2(:,e)'*(r2) + q3(:,e)'*(r3) + q4(:,e)'*(r4)));
    %         UavgWADG(i) = UavgWADG(i) + sum(M*((Vq'*diag(wq./J(:,e))*Vq)\(M*Pq*(rho(:,e)+rhou(:,e)+rhov(:,e)+E(:,e)))));
    %         rhsavgWADG(i) = rhsavgWADG(i) + sum(M*((Vq'*diag(wq./J(:,e))*Vq)\(M*Pq*(rhs1(:,e)+rhs2(:,e)+rhs3(:,e)+rhs4(:,e)))));
    %     end
    
%         rhsavg(i) = abs(sum(sum(wJq.*(rhs1)))) + abs(sum(sum(wJq.*(rhs2)))) + abs(sum(sum(wJq.*(rhs3)))) + abs(sum(sum(wJq.*(rhs4))));
%         rhsavg(i) = abs(sum(sum(wJq.*(rhs1+rhs2+rhs3+rhs4))));
    %     Uavg(i) = sum(sum(wJq.*(rho+rhou+rhov+E)));    
    
%     Sq = -rho.*s(rho,rhou,rhov,E);
%     entropy(i) = sum(sum(wJq.*Sq));

if  (mod(i,10)==0 || i==Nsteps)
    clf
    %     pp = rho;
    %     vv = real(Vp*Pq*pp);
    %     color_line3(xp,yp,vv,vv,'.');
    PlotField2D(15,x,y,-Pq*rho);view(2)
    axis equal
    axis tight
    colorbar
    title(sprintf('time = %f, N = %d, K1D = %d. Step %d/%d',dt*i,N,K1D,i,Nsteps))
    %title(sprintf('time = %f, N = %d, K1D = %d, entropy rhs = %g',dt*i,N,K1D,rhstest(i)))
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

return

figure(2)
% plot((1:Nsteps)*dt,abs(UavgWADG-UavgWADG(1)),'o--')
hold on
% plot((1:Nsteps)*dt,abs(Uavg-Uavg(1)),'x--','DisplayName','Solution average')
% plot((1:Nsteps)*dt,rhsavgWADG,'s--','DisplayName','RHS WADG avg')
plot((1:Nsteps)*dt,rhsavg,'^--','DisplayName','RHS avg')
legend show
%legend('WADG average','Average','RHS WADG avg','RHS avg')
xlabel('Time')
ylabel('Conservation error')
set(gca,'fontsize',16)

% dn = round(Nsteps/25);
% print_pgf_coordinates((1:dn:Nsteps)*dt,rhsavg(1:dn:Nsteps))


function [rhs1 rhs2 rhs3 rhs4] = RHS2Dsimple(rhoq,uq,vq,Eq,rhoM,uM,vM,EM,QM,time)

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

% Lax-Friedrichs flux
global gamma tau
unorm2 = (uM.^2+vM.^2);
pM = (gamma-1)*(EM - .5*rhoM.*unorm2);

% % % walls = top/bottom
% global mapBwall
% uP(mapBwall) = uM(mapBwall);
% vP(mapBwall) = -vM(mapBwall);
% rhoP(mapBwall) = rhoM(mapBwall);
% pP = pM(mapPq);
% pP(mapBwall) = pM(mapBwall);
% EP(mapBwall) = pP(mapBwall)/(gamma-1) + .5*rhoP(mapBwall).*(uP(mapBwall).^2+vP(mapBwall).^2); 

QM{1} = rhoM;
QM{2} = rhoM.*uM;
QM{3} = rhoM.*vM;
QM{4} = EM;

QP{1} = rhoP;
QP{2} = rhoP.*uP;
QP{3} = rhoP.*vP;
QP{4} = EP;

% QP{1} = QM{1}(mapPq); QP{2} = QM{2}(mapPq);
% QP{3} = QM{3}(mapPq); QP{4} = QM{4}(mapPq);

if 1
    cvel = sqrt(gamma*pM./rhoM);
    lam = sqrt(unorm2)+cvel;
    LFc = max(lam(mapPq),lam);

    dQ1 = QP{1}-QM{1};
    dQ2 = QP{2}-QM{2};
    dQ3 = QP{3}-QM{3};
    dQ4 = QP{4}-QM{4};
    Lf1 = tau*LFc.*dQ1.*sJ;
    Lf2 = tau*LFc.*dQ2.*sJ;
    Lf3 = tau*LFc.*dQ3.*sJ;
    Lf4 = tau*LFc.*dQ4.*sJ;
else
    global V1 V2 V3 V4
    du1 = (V1(QP{1},QP{2},QP{3},QP{4})-V1(QM{1},QM{2},QM{3},QM{4}));
    du2 = (V2(QP{1},QP{2},QP{3},QP{4})-V2(QM{1},QM{2},QM{3},QM{4}));
    du3 = (V3(QP{1},QP{2},QP{3},QP{4})-V3(QM{1},QM{2},QM{3},QM{4}));
    du4 = (V4(QP{1},QP{2},QP{3},QP{4})-V4(QM{1},QM{2},QM{3},QM{4}));
            
    u2avg = avg(uM.^2,uP.^2);
    v2avg = avg(vM.^2,vP.^2);
    uavg = avg(uM,uP);
    vavg = avg(vM,vP);
    betaM = beta(rhoM,uM,vM,EM);
    betaP = beta(rhoP,uP,vP,EP);
    betalog = logmean(betaM,betaP);
    pavg = avg(rhoM,rhoP)./(2*avg(betaM,betaP));   
    rholog = logmean(rhoM,rhoP); 
    
    ubar = 2*(uavg.^2+vavg.^2)-(u2avg+v2avg);
    un = (uavg.*nx+vavg.*ny);
    h = gamma./(2*betalog*(gamma-1)) + .5*ubar;
    a = sqrt(gamma.*pavg./rholog);
    
    R11 = 1; R12 = 1; R13 = 0; R14 = 1;
    R21 = uavg-a.*nx; R22 = uavg; R23 = ny; R24  = uavg+a.*nx;
    R31 = vavg-a.*ny; R32 = vavg; R33 = -nx; R34  = vavg+a.*ny;
    R41 = h-a.*un; R42 = .5*ubar; R43 = uavg.*ny - vavg.*nx; R44 = h+a.*un;
    D11 = abs((un - a)).*rholog/(2*gamma);
    D22 = abs((un)).*rholog*(gamma-1)/gamma;
    D33 = abs((un)).*pavg;
    D44 = abs((un + a)).*rholog/(2*gamma);
    
%         R11 = 1; R12 = 1; R13 = 0; R14 = 1;
%         R21 = uavg-a; R22 = uavg; R23 = 0; R24  = uavg+a;
%         R31 = vavg; R32 = vavg; R33 = 1; R34  = vavg;
%         R41 = h-a.*uavg; R42 = .5*ubar; R43 = vavg; R44 = h+a.*uavg;
%         D11 = abs((uavg - a)).*rholog/(2*gamma);
%         D22 = abs((uavg)).*rholog*(gamma-1)/gamma;
%         D33 = abs((uavg)).*pavg;
%         D44 = abs((uavg + a)).*rholog/(2*gamma);

    r1 = D11.*(R11.*du1 + R21.*du2 + R31.*du3 + R41.*du4);
    r2 = D22.*(R12.*du1 + R22.*du2 + R32.*du3 + R42.*du4);
    r3 = D33.*(R13.*du1 + R23.*du2 + R33.*du3 + R43.*du4);
    r4 = D44.*(R14.*du1 + R24.*du2 + R34.*du3 + R44.*du4);
    
    Lf1 = tau*(R11.*r1 + R12.*r2 + R13.*r3 + R14.*r4).*sJ;
    Lf2 = tau*(R21.*r1 + R22.*r2 + R23.*r3 + R24.*r4).*sJ;
    Lf3 = tau*(R31.*r1 + R32.*r2 + R33.*r3 + R34.*r4).*sJ;
    Lf4 = tau*(R41.*r1 + R42.*r2 + R43.*r3 + R44.*r4).*sJ;
end

fSf1 = nxJ.*fxS1(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + nyJ.*fyS1(rhoM,uM,vM,EM,rhoP,uP,vP,EP);
fSf2 = nxJ.*fxS2(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + nyJ.*fyS2(rhoM,uM,vM,EM,rhoP,uP,vP,EP);
fSf3 = nxJ.*fxS3(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + nyJ.*fyS3(rhoM,uM,vM,EM,rhoP,uP,vP,EP);
fSf4 = nxJ.*fxS4(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + nyJ.*fyS4(rhoM,uM,vM,EM,rhoP,uP,vP,EP);

fSf1 = fSf1  - .5*Lf1;
fSf2 = fSf2  - .5*Lf2;
fSf3 = fSf3  - .5*Lf3;
fSf4 = fSf4  - .5*Lf4;

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
rxJK = rxJ(1,e);
sxJK = sxJ(1,e);
ryJK = ryJ(1,e);
syJK = syJ(1,e);
    
    Dx = DNr.*rxJK + DNs.*sxJK;
    Dy = DNr.*ryJK + DNs.*syJK;
    
    divF1(:,e) = sum(Dx.*FxS1,2) + sum(Dy.*FyS1,2);
    divF2(:,e) = sum(Dx.*FxS2,2) + sum(Dy.*FyS2,2);
    divF3(:,e) = sum(Dx.*FxS3,2) + sum(Dy.*FyS3,2);
    divF4(:,e) = sum(Dx.*FxS4,2) + sum(Dy.*FyS4,2);

end

rhs1 = 2*[VqPq VqLq]*divF1 + VqLq*(fSf1);
rhs2 = 2*[VqPq VqLq]*divF2 + VqLq*(fSf2);
rhs3 = 2*[VqPq VqLq]*divF3 + VqLq*(fSf3);
rhs4 = 2*[VqPq VqLq]*divF4 + VqLq*(fSf4);

% apply wadg and correction
% global wq
% a1 = (wq'*rhs1); % store averages
% a2 = (wq'*rhs2);
% a3 = (wq'*rhs3);
% a4 = (wq'*rhs4);
% rhs1 = VqPq*(rhs1./J);
% rhs2 = VqPq*(rhs2./J);
% rhs3 = VqPq*(rhs3./J);
% rhs4 = VqPq*(rhs4./J);
% rhs1 = -rhs1;
% rhs2 = -rhs2;
% rhs3 = -rhs3;
% rhs4 = -rhs4;

rhs1 = -(rhs1./J);
rhs2 = -(rhs2./J);
rhs3 = -(rhs3./J);
rhs4 = -(rhs4./J);

end

function [rho u v p] = vortexSolution(x,y,t)

global gamma
xs = .5;
x0 = 0+.25;
y0 = .5;

% shock
Ms = 1.1;

rho_up = 1;
u_up = Ms*sqrt(gamma);
v_up = 0;
p_up = 1;

rho_ratio_ud = (2+(gamma-1)*Ms^2)/((gamma+1)*Ms^2);
u_ratio_ud = 1/rho_ratio_ud;
p_ratio_ud = 1/(1 + 2*gamma/(gamma+1)*(Ms^2-1));

% uup = mean(x)<xs;
rhos = rho_up + 0*x;
us = u_up + 0*x;
ps = p_up + 0*x;
vs = 0*x;

ds = mean(x)>xs;
rhos(:,ds) = rho_up / rho_ratio_ud;
us(:,ds) = u_up / u_ratio_ud;
ps(:,ds) = p_up / p_ratio_ud;

% % vortex
% Mv = .9;
% vm = Mv*sqrt(gamma);
% a = .075; b = .175;
% vtheta = (vm * r / a).*(r<=a) + (vm * (a/(a^2-b^2)).*(r-b^2./r)).*(r > a).*(r<b);

% from shu's paper "Efficient Implementation of Weighted ENO Schemes"
r = sqrt((x-x0).^2+(y-y0).^2);
epsilon = .3; 
alpha = .204;
rc = .05; 
tau = r./rc;
vtheta = epsilon * tau .* exp(alpha*(1-tau.^2));
theta = atan2(y-y0,x-x0);

Tu = p_up/rho_up;
Tvor = (p_up/rho_up) - (gamma-1)*epsilon^2.*exp(2*alpha*(1-tau.^2))/(4*alpha*gamma);

% u = us;
% v = vs;
% rho = rhos;
% p = ps;

u = us + vtheta.*sin(theta);
v = vs - vtheta.*cos(theta);
rho = rhos .* (Tvor/Tu).^((1/(gamma-1)));
p = ps .* (Tvor/Tu).^((1/(gamma-1)));


end