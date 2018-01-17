clear
Globals2D

N = 3;
K1D = 8;
FinalTime = 2.5;
CFL = .125;
global tau
tau = 1;
projectV = 0;
a = 1/8;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
% VX = VX/max(abs(VX));  VY = VY/max(abs(VY));
% VX = VX*5 + 5; VY = VY*5;

% iids = abs(VX) < 1 & abs(VY) < 1;
% VX(iids) = VX(iids) + .25/K1D*randn(size(VX(iids)));
% VY(iids) = VY(iids) + .25/K1D*randn(size(VY(iids)));

StartUp2D;
BuildPeriodicMaps2D(max(VX)-min(VX),max(VY)-min(VY));
% StartUp2D;

% plotting nodes
[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
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

x = x + a*cos(pi/2*x).*cos(3*pi/2*y);
y = y + a*cos(3*pi/2*x).*cos(pi/2*y);

xq = Vq*x; yq = Vq*y;
xp = Vp*x; yp = Vp*y;
plot(x(Fmask(:),:),y(Fmask(:),:),'k-','linewidth',2);return

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

x0 = mean(VX); y0 = 0;

hex = @(x,y) 2 + exp(-50*((x-x0).^2 + (y-y0).^2));
% hex = @(x,y) 2 + exp(-5^2*((x-x0)-(y-y0)).^2);
% hex = @(x,y) 2 + exp(-1^2*(x-x0).^2);
% hex = @(x,y) 2 + (abs(x)<.5).*(abs(y)<.5);
h = Pq*hex(xq,yq);
hu = Pq*(xq*0);
hv = Pq*(xq*0);

% bathymetry + gravity
b = 0;
g = 1;

global fxS1 fxS2 fxS3 fyS1 fyS2 fyS3
global avg
avg = @(x,y) .5*(x+y);
if projectV
    fxS1 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(uL,uR);
    fxS2 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(uL,uR).^2 + .5*g*avg(hL.^2,hR.^2);
    fxS3 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(uL,uR).*avg(vL,vR);
    
    fyS1 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(vL,vR);
    fyS2 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(uL,uR).*avg(vL,vR);
    fyS3 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(vL,vR).^2 + .5*g*avg(hL.^2,hR.^2);
else
    fxS1 = @(hL,uL,vL,hR,uR,vR) avg(hL.*uL,hR.*uR);
    fxS2 = @(hL,uL,vL,hR,uR,vR) avg(hL.*uL.^2,hR.*uR.^2) + .5*g*avg(hL.^2,hR.^2);
    fxS3 = @(hL,uL,vL,hR,uR,vR) avg(hL.*uL.*vL,hR.*uR.*vR);
    
    fyS1 = @(hL,uL,vL,hR,uR,vR) avg(hL.*vL,hR.*vR);
    fyS2 = @(hL,uL,vL,hR,uR,vR) avg(hL.*uL.*vL,hR.*uR.*vR);
    fyS3 = @(hL,uL,vL,hR,uR,vR) avg(hL.*vL.^2,hR.*vR.^2) + .5*g*avg(hL.^2,hR.^2);
end

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
        
        if projectV
            % project to entropy variables
            hq = Vq*h;
            huq = Vq*hu;
            hvq = Vq*hv;
            uq = huq./hq;
            vq = hvq./hq;
            
            q1 = Pq*(g.*(hq+b)- .5*(uq.^2 + vq.^2));
            u  = Pq*uq;
            v  = Pq*vq;
            
            % evaluate at quad/surface points
            qq = Vq*q1; qM = Vfq*q1;
            uq = Vq*u;  uM = Vfq*u;
            vq = Vq*v;  vM = Vfq*v;
            hq = (qq + .5*(uq.^2 + vq.^2))./g - b;
            hM = (qM + .5*(uM.^2 + vM.^2))./g - b;
        else
            hq = Vq*h;
            huq = Vq*hu;
            hvq = Vq*hv;
            uq = huq./hq;
            vq = hvq./hq;
            
            hM = Vfq*h;
            uM = (Vfq*hu)./hM;
            vM = (Vfq*hv)./hM;
            qq = g.*(hq+b)- .5*(uq.^2 + vq.^2);
            uq = huq./hq;
            vq = hvq./hq;
        end
                
        [rhs1 rhs2 rhs3]  = RHS2D(hq,uq,vq,hM,uM,vM);
        
        if INTRK==5
            rhstest(i) = sum(sum(wJq.*((Vq*rhs1).*qq + (Vq*rhs2).*uq + (Vq*rhs3).*vq)));
        end
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        
        h  = h  + rk4b(INTRK)*res1;
        hu = hu + rk4b(INTRK)*res2;
        hv = hv + rk4b(INTRK)*res3;                
    end;
    
    % project to entropy variables
    hq = Vq*h;
    huq = Vq*hu;
    hvq = Vq*hv;
    uq = huq./hq;
    vq = hvq./hq;
    energyq = (.5*hq.*(uq.^2+vq.^2) + .5*g*hq.^2 + g.*hq.*b);
    energy(i) = sum(sum(wJq.*energyq));
    
    if mod(i,5)==0 || i==Nsteps
        clf
        pp = h;
        vv = real(Vp*pp);
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f',dt*i))
                view(3)
        drawnow
        
    end
    
end

% return
figure(2)
semilogy(dt*(1:i),energy,'--','linewidth',2)
% dS = abs(energy-energy(1));
% semilogy(dt*(1:Nsteps),dS,'--','linewidth',2)
hold on
% semilogy(dt*(1:Nsteps),abs(rhstest),'x')


function [rhs1 rhs2 rhs3] = RHS2D(hq,uq,vq,hM,uM,vM)

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

hP = hM(mapPq);
uP = uM(mapPq);
vP = vM(mapPq);

% Lax-Friedrichs flux
g = 1;
cvel = sqrt(g*hM);
lfm = sqrt(uM.^2+vM.^2)+cvel;
LFc = max(lfm(mapPq),lfm);
hUn = (hP.*uP-hM.*uM).*nx + (hP.*vP-hM.*vM).*ny;
global tau
Lf1 = tau*LFc.*(hP-hM).*sJ;
Lf2 = tau*LFc.*(hUn.*nx).*sJ;
Lf3 = tau*LFc.*(hUn.*ny).*sJ;

fSf1 = nxJ.*fxS1(hM,uM,vM,hP,uP,vP) + nyJ.*fyS1(hM,uM,vM,hP,uP,vP);
fSf2 = nxJ.*fxS2(hM,uM,vM,hP,uP,vP) + nyJ.*fyS2(hM,uM,vM,hP,uP,vP);
fSf3 = nxJ.*fxS3(hM,uM,vM,hP,uP,vP) + nyJ.*fyS3(hM,uM,vM,hP,uP,vP);

rhs1 = zeros(Np,K);
rhs2 = zeros(Np,K);
rhs3 = zeros(Np,K);
for e = 1:K
    
    % local arithmetic operations - form on the fly for GPU
    [hx hy] = meshgrid(hq(:,e));  [hfx hfy] = meshgrid(hM(:,e),hq(:,e));
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

