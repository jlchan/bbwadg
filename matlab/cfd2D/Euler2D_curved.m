clear
Globals2D

N = 3;
K1D = 8;
FinalTime = 1;
CFL = .75;
global tau
tau = 1;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

StartUp2D;
BuildPeriodicMaps2D(2,2);
% StartUp2D;

% plotting nodes
[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;
% plot(xp,yp,'o')
% return

global Vq Pq Lq Lqf Vfqf Vfq Pfqf
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ nrJq nsJq
global mapPq
Nq = 2*N;
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
                %                     keyboard
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

a = 1/8;
x = x + a*sin(pi*x).*sin(pi*y);
y = y + a*sin(pi*x).*sin(pi*y);

xq = Vq*x; yq = Vq*y;
xp = Vp*x; yp = Vp*y;
% plot(x(Fmask(:),:),y(Fmask(:),:),'.');return

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


%% problem params setup

x0 = 0; y0 = 0;

hex = @(x,y) 2 + exp(-5^2*(x-x0).^2);
rho = Pq*hex(xq,yq);
rhou = Pq*(xq*0);
rhov = Pq*(xq*0);

% bathymetry + gravity
g = 1;

rhoe = @(rho,rhou,rhov,E) E - (rhou.^2+rhov.^2)./(2*rho);
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


global fxS1 fxS2 fxS3 fyS1 fyS2 fyS3
global avg
avg = @(x,y) .5*(x+y);
fxS1 = @(rhoL,uL,vL,rhoR,uR,vR) logmean(rhoL,rhoR).*avg(uL,uR);
fxS2 = @(rhoL,uL,vL,rhoR,uR,vR) avg(rhoL,rhoR).*avg(uL,uR).^2 + .5*g*avg(rhoL.^2,rhoR.^2);
fxS3 = @(rhoL,uL,vL,rhoR,uR,vR) avg(rhoL,rhoR).*avg(uL,uR).*avg(vL,vR);

fyS1 = @(rhoL,uL,vL,rhoR,uR,vR) avg(rhoL,rhoR).*avg(vL,vR);
fyS2 = @(rhoL,uL,vL,rhoR,uR,vR) avg(rhoL,rhoR).*avg(uL,uR).*avg(vL,vR);
fyS3 = @(rhoL,uL,vL,rhoR,uR,vR) avg(rhoL,rhoR).*avg(vL,vR).^2 + .5*g*avg(rhoL.^2,rhoR.^2);




%%

global wJq
wJq = diag(wq)*(J);

% Runge-Kutta residual storage
res1 = zeros(Np,K); res2 = zeros(Np,K); res3 = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
CNh = max(CN*max(sJ(:)./Jf(:)));
dt = CFL*2/CNh;

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

figure(1)
for i = 1:Nsteps
    for INTRK = 1:5
        
        % project to entropy variables        
        rhoq = Vq*rho;
        rhouq = Vq*rhou;
        rhovq = Vq*rhov;
        uq = rhouq./rhoq;
        vq = rhovq./rhoq;
        
        q1 = Pq*(g.*(rhoq+b)- .5*(uq.^2 + vq.^2));
        u  = Pq*uq;
        v  = Pq*vq;
        
        % evaluate at quad/surface points
        qq = Vq*q1; qM = Vfq*q1;
        uq = Vq*u;  uM = Vfq*u;
        vq = Vq*v;  vM = Vfq*v;
        rhoq = (qq + .5*(uq.^2 + vq.^2))./g - b;
        rhoM = (qM + .5*(uM.^2 + vM.^2))./g - b;
                
        [rhs1 rhs2 rhs3]  = RHS2D(rhoq,uq,vq,rhoM,uM,vM);
        
        if INTRK==5
            rhstest(i) = sum(sum(wJq.*((Vq*rhs1).*qq + (Vq*rhs2).*uq + (Vq*rhs3).*vq)));
        end
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        
        rho  = rho  + rk4b(INTRK)*res1;
        rhou = rhou + rk4b(INTRK)*res2;
        rhov = rhov + rk4b(INTRK)*res3;                
    end;
    
    % project to entropy variables
    rhoq = Vq*rho;
    rhouq = Vq*rhou;
    rhovq = Vq*rhov;
    uq = rhouq./rhoq;
    vq = rhovq./rhoq;
    energy(i) = sum(sum(wJq.*(.5*rhoq.*(uq.^2+vq.^2) + .5*g*rhoq.^2 + g.*rhoq.*b)));
    
    if mod(i,5)==0 || i==Nsteps
        clf
        pp = rho;
        vv = Vp*pp;
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f',dt*i))
        %         view(3)
        drawnow
        
    end
    
end

% return
figure(2)
semilogy(dt*(1:Nsteps),energy,'--')
hold on
semilogy(dt*(1:Nsteps),abs(rhstest),'x')


function [rhs1 rhs2 rhs3] = RHS2D(rhoq,uq,vq,rhoM,uM,vM)

Globals2D;

global Vq Pq Lq Lqf Vfqf Vfq Pfqf
% global Drq Dsq Lrq Lsq
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ nrJq nsJq
global mapPq
global fxS1 fxS2 fxS3 fyS1 fyS2 fyS3
global avg

% operators
Drq = (Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vfq*Pq);
Dsq = (Vq*Ds*Pq - .5*Vq*Lq*diag(nsJ)*Vfq*Pq);
VfPq = (Vfq*Pq);
VqLq = Vq*Lq;

rhoP = rhoM(mapPq);
uP = uM(mapPq);
vP = vM(mapPq);

% Lax-Friedrichs flux
g = 1;
cvel = sqrt(g*rhoM);
LFc = max(sqrt(uM.^2+vM.^2)+cvel);
rhoun = (rhoP.*uP-rhoM.*uM).*nx + (rhoP.*vP-rhoM.*vM).*ny;
global tau
Lf1 = tau*LFc.*(rhoP-rhoM).*sJ;
Lf2 = tau*LFc.*(rhoun.*nx).*sJ;
Lf3 = tau*LFc.*(rhoun.*ny).*sJ;

fSf1 = nxJ.*fxS1(rhoM,uM,vM,rhoP,uP,vP) + nyJ.*fyS1(rhoM,uM,vM,rhoP,uP,vP);
fSf2 = nxJ.*fxS2(rhoM,uM,vM,rhoP,uP,vP) + nyJ.*fyS2(rhoM,uM,vM,rhoP,uP,vP);
fSf3 = nxJ.*fxS3(rhoM,uM,vM,rhoP,uP,vP) + nyJ.*fyS3(rhoM,uM,vM,rhoP,uP,vP);

rhs1 = zeros(Np,K);
rhs2 = zeros(Np,K);
rhs3 = zeros(Np,K);
for e = 1:K
    
    % local aritrhoMetic operations - form on the fly for GPU
    [hx hy] = meshgrid(rhoq(:,e));  [hfx hfy] = meshgrid(rhoM(:,e),rhoq(:,e));
    [ux uy] = meshgrid(uq(:,e));  [ufx ufy] = meshgrid(uM(:,e),uq(:,e));
    [vx vy] = meshgrid(vq(:,e));  [vfx vfy] = meshgrid(vM(:,e),vq(:,e));
    
    % avoiding geometric aliasing
    [rxJ1 rxJ2] = meshgrid(rxJ(:,e));  [sxJ1 sxJ2] = meshgrid(sxJ(:,e));
    [ryJ1 ryJ2] = meshgrid(ryJ(:,e));  [syJ1 syJ2] = meshgrid(syJ(:,e));
    rxJK = avg(rxJ1,rxJ2);  sxJK = avg(sxJ1,sxJ2);
    ryJK = avg(ryJ1,ryJ2);  syJK = avg(syJ1,syJ2);
    
    FxS1 = fxS1(hx,ux,vx,hy,uy,vy);  FyS1 = fyS1(hx,ux,vx,hy,uy,vy);      
    FxS2 = fxS2(hx,ux,vx,hy,uy,vy);  FyS2 = fyS2(hx,ux,vx,hy,uy,vy);  
    FxS3 = fxS3(hx,ux,vx,hy,uy,vy);  FyS3 = fyS3(hx,ux,vx,hy,uy,vy);  
        
    % premultiply by geofacs
    FrS1 = rxJK.*FxS1 + ryJK.*FyS1; 
    FsS1 = sxJK.*FxS1 + syJK.*FyS1;
    FrS2 = rxJK.*FxS2 + ryJK.*FyS2; 
    FsS2 = sxJK.*FxS2 + syJK.*FyS2;
    FrS3 = rxJK.*FxS3 + ryJK.*FyS3; 
    FsS3 = sxJK.*FxS3 + syJK.*FyS3;        
        
    % bulk of GPU work: application of local operators
    divF1 = sum(Drq.*FrS1,2) + sum(Dsq.*FsS1,2);
    divF2 = sum(Drq.*FrS2,2) + sum(Dsq.*FsS2,2);
    divF3 = sum(Drq.*FrS3,2) + sum(Dsq.*FsS3,2);        
    
    %% flux/volume combos
    
    FxSf1 = fxS1(hfx,ufx,vfx,hfy,ufy,vfy);  FySf1 = fyS1(hfx,ufx,vfx,hfy,ufy,vfy);
    FxSf2 = fxS2(hfx,ufx,vfx,hfy,ufy,vfy);  FySf2 = fyS2(hfx,ufx,vfx,hfy,ufy,vfy);
    FxSf3 = fxS3(hfx,ufx,vfx,hfy,ufy,vfy);  FySf3 = fyS3(hfx,ufx,vfx,hfy,ufy,vfy);
    
    % averaging to prevent geometric aliasing?     
    [rxJ1 rxJ2] = meshgrid(rxJ(:,e),rxJf(:,e)); rxJK = avg(rxJ1,rxJ2)';
    [sxJ1 sxJ2] = meshgrid(sxJ(:,e),sxJf(:,e)); sxJK = avg(sxJ1,sxJ2)';
    [ryJ1 ryJ2] = meshgrid(ryJ(:,e),ryJf(:,e)); ryJK = avg(ryJ1,ryJ2)';
    [syJ1 syJ2] = meshgrid(syJ(:,e),syJf(:,e)); syJK = avg(syJ1,syJ2)';
    
    Nq = size(FxSf1,1);
    FxSf1r = FxSf1.*nrJq; FxSf1s = FxSf1.*nsJq;         
    FxSf2r = FxSf2.*nrJq; FxSf2s = FxSf2.*nsJq;  
    FxSf3s = FxSf3.*nsJq; FxSf3r = FxSf3.*nrJq; 
    
    FySf1r = FySf1.*nrJq; FySf1s = FySf1.*nsJq;        
    FySf2r = FySf2.*nrJq; FySf2s = FySf2.*nsJq;    
    FySf3s = FySf3.*nsJq; FySf3r = FySf3.*nrJq;
    
    FSf1 = rxJK.*FxSf1r + ryJK.*FySf1r + sxJK.*FxSf1s + syJK.*FySf1s;
    FSf2 = rxJK.*FxSf2r + ryJK.*FySf2r + sxJK.*FxSf2s + syJK.*FySf2s;
    FSf3 = rxJK.*FxSf3r + ryJK.*FySf3r + sxJK.*FxSf3s + syJK.*FySf3s;
        
    % bulk of GPU work: application of local operators
    fSproj1 = sum(VfPq.*FSf1',2);
    fSproj2 = sum(VfPq.*FSf2',2);
    fSproj3 = sum(VfPq.*FSf3',2);     
    
    % intermediate fluxes - form on the fly on GPU
    f1 = FSf1 + repmat((fSf1(:,e) - fSproj1 - .25*Lf1(:,e))',Nq,1) ;
    f2 = FSf2 + repmat((fSf2(:,e) - fSproj2 - .25*Lf2(:,e))',Nq,1) ;
    f3 = FSf3 + repmat((fSf3(:,e) - fSproj3 - .25*Lf3(:,e))',Nq,1) ;
    
    %% project back to polynomial space - can incorporate WADG here
    
    rhs1(:,e) =  Pq*(divF1 + .5*sum(VqLq.*f1,2));
    rhs2(:,e) =  Pq*(divF2 + .5*sum(VqLq.*f2,2));
    rhs3(:,e) =  Pq*(divF3 + .5*sum(VqLq.*f3,2));

end

rhs1 = -2*Pq*((Vq*rhs1)./J);
rhs2 = -2*Pq*((Vq*rhs2)./J);
rhs3 = -2*Pq*((Vq*rhs3)./J);

end

