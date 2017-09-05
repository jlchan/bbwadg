clear
Globals2D

N = 4;
K1D = 16;
FinalTime = 1.5;
CFL = .75;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
% iids = abs(VX) < 1 & abs(VY) < 1;
% VX(iids) = VX(iids) + .25/K1D*randn(size(VX(iids)));
% VY(iids) = VY(iids) + .25/K1D*randn(size(VY(iids)));

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

if 0
    u = (1:Np)';
    uq = Vq*u;
    uf = Vfq*u;
    
    [ux uy] = meshgrid(uq);
    [ufx ufy] = meshgrid(uf,uq);
    
    FS = fS(ux,uy);
    FSf = fS(ufx,ufy);
    
    Pq*sum((diag(wq)*Vq*Lq).*FSf,2)
    sum((diag(wq)*Vq*Lq).*FSf,1)*(M\Vfq')'
    M\Vfq'*sum((diag(wfq)*Vfq*Pq).*FSf',2)
end

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
                
                %                     plot(x1,y1,'o')
                %                     hold on
                %                     plot(x2,y2,'x')
                %                     axis([-1,1,-1,1])
                %                     keyboard
                
            end
            mapPq(id1,e) = id2(p) + (enbr-1)*(Nfq*Nfaces);
        end
    end
end

% hold on
% for i = 1:nnz(xf)
%     id = mapPq(i) ;
%     plot(xf(i),yf(i),'o')
%     plot(xf(id),yf(id),'x')
%    pause
% end
% return

%% make curvilinear mesh (still unstable?)

a = 1/16;
x = x + a*sin(pi*x).*sin(pi*y);
y = y + a*sin(pi*x).*sin(pi*y);

xq = Vq*x; yq = Vq*y;
xp = Vp*x; yp = Vp*y;
% plot(xp,yp,'.');return

rxJ = zeros(Nq,K); sxJ = zeros(Nq,K);
ryJ = zeros(Nq,K); syJ = zeros(Nq,K);
J = zeros(Nq,K);
rxJf = zeros(Nfq*Nfaces,K); sxJf = zeros(Nfq*Nfaces,K);
ryJf = zeros(Nfq*Nfaces,K); syJf = zeros(Nfq*Nfaces,K);
for e = 1:K
    [rxk,sxk,ryk,syk,Jk] = GeometricFactors2D(x(:,e),y(:,e),Vq*Dr,Vq*Ds);
    rxJ(:,e) = rxk.*Jk;    sxJ(:,e) = sxk.*Jk;
    ryJ(:,e) = ryk.*Jk;    syJ(:,e) = syk.*Jk;
    J(:,e) = Jk;
    
    [rxk,sxk,ryk,syk,Jfk] = GeometricFactors2D(x(:,e),y(:,e),Vfq*Dr,Vfq*Ds);
    rxJf(:,e) = rxk.*Jfk;    sxJf(:,e) = sxk.*Jfk;
    ryJf(:,e) = ryk.*Jfk;    syJf(:,e) = syk.*Jfk;
end

%% problem params setup

x0 = 0; y0 = 0;

%hex = @(x,y) 2+exp(-5^2*((x-x0).^2 + (y-y0).^2));
hex = @(x,y) 2 + exp(-5^2*(x-x0).^2);
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
fxS1 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(uL,uR);
fxS2 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(uL,uR).^2 + .5*g*avg(hL.^2,hR.^2);
fxS3 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(uL,uR).*avg(vL,vR);

fyS1 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(vL,vR);
fyS2 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(uL,uR).*avg(vL,vR);
fyS3 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(vL,vR).^2 + .5*g*avg(hL.^2,hR.^2);

%%

global wJq
wJq = diag(wq)*(J);

% Runge-Kutta residual storage
res1 = zeros(Np,K); res2 = zeros(Np,K); res3 = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
CNh = max(CN*max(Fscale(:)));
dt = CFL*2/CNh;

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

figure(1)
for i = 1:Nsteps
    for INTRK = 1:5
        
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
                
        [rhs1 rhs2 rhs3]  = RHS2D(hq,uq,vq,hM,uM,vM);
        
        if INTRK==5
            rhstest(i) = sum(sum(wJq.*((Vq*rhs1).*qq + (Vq*rhs2).*uq + (Vq*rhs3).*vq)));
        end
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        
        h = h + rk4b(INTRK)*res1;
        hu = hu + rk4b(INTRK)*res2;
        hv = hv + rk4b(INTRK)*res3;                
    end;
    
    %     uq = Vq*u;
    %     energy(i) = sum(sum(wJq.*uq.^2));
    
    if mod(i,5)==0 || i==Nsteps
        clf
        pp = h;
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

return
figure(2)
semilogy(dt*(1:Nsteps),energy,'--')
hold on
semilogy(dt*(1:Nsteps),abs(rhstest),'x')


function [rhs1 rhs2 rhs3] = RHS2D(hq,uq,vq,hM,uM,vM)

Globals2D;

global Vq Pq Lq Lqf Vfqf Vfq Pfqf
% global Drq Dsq Lrq Lsq
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ nrJq nsJq
global mapPq
global fxS1 fxS2 fxS3 fyS1 fyS2 fyS3
global avg

% quadrature-based operators 
Drq = (Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vfq*Pq);
Dsq = (Vq*Ds*Pq - .5*Vq*Lq*diag(nsJ)*Vfq*Pq);
VfPq = (Vfq*Pq);
VqLq = Vq*Lq;
Nq = size(Vq,1);

% extract + values
hP = hM(mapPq);
uP = uM(mapPq);
vP = vM(mapPq);

% precompute normal fluxes
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
        
    rxJK = rxJ(:,e); sxJK = sxJ(:,e);
    ryJK = ryJ(:,e); syJK = syJ(:,e);    
    
    FxS1 = fxS1(hx,ux,vx,hy,uy,vy);  FyS1 = fyS1(hx,ux,vx,hy,uy,vy);      
    FxS2 = fxS2(hx,ux,vx,hy,uy,vy);  FyS2 = fyS2(hx,ux,vx,hy,uy,vy);  
    FxS3 = fxS3(hx,ux,vx,hy,uy,vy);  FyS3 = fyS3(hx,ux,vx,hy,uy,vy);  
        
    % premultiply by geofacs
%     FrS1 = rxJK.*FxS1 + ryJK.*FyS1;  FsS1 = sxJK.*FxS1 + syJK.*FyS1;
%     FrS2 = rxJK.*FxS2 + ryJK.*FyS2;  FsS2 = sxJK.*FxS2 + syJK.*FyS2;
%     FrS3 = rxJK.*FxS3 + ryJK.*FyS3;  FsS3 = sxJK.*FxS3 + syJK.*FyS3;        
        
    % bulk of GPU work: application of local operators
    divF1 = rxJK.*sum(Drq.*FxS1,2) + sxJK.*sum(Dsq.*FxS1,2) + ryJK.*sum(Drq.*FyS1,2) + syJK.*sum(Dsq.*FyS1,2);
    divF2 = rxJK.*sum(Drq.*FxS2,2) + sxJK.*sum(Dsq.*FxS2,2) + ryJK.*sum(Drq.*FyS2,2) + syJK.*sum(Dsq.*FyS2,2);
    divF3 = rxJK.*sum(Drq.*FxS3,2) + sxJK.*sum(Dsq.*FxS3,2) + ryJK.*sum(Drq.*FyS3,2) + syJK.*sum(Dsq.*FyS3,2);
%     divF2 = sum(Drq.*FrS2,2) + sum(Dsq.*FsS2,2);
%     divF3 = sum(Drq.*FrS3,2) + sum(Dsq.*FsS3,2);
    
    %% flux/volume combos
    
    FxSf1 = fxS1(hfx,ufx,vfx,hfy,ufy,vfy);  FySf1 = fyS1(hfx,ufx,vfx,hfy,ufy,vfy);
    FxSf2 = fxS2(hfx,ufx,vfx,hfy,ufy,vfy);  FySf2 = fyS2(hfx,ufx,vfx,hfy,ufy,vfy);
    FxSf3 = fxS3(hfx,ufx,vfx,hfy,ufy,vfy);  FySf3 = fyS3(hfx,ufx,vfx,hfy,ufy,vfy);
    
    FSf1 = FxSf1.*repmat(nxJ(:,e)',Nq,1) + FySf1.*repmat(nyJ(:,e)',Nq,1);
    FSf2 = FxSf2.*repmat(nxJ(:,e)',Nq,1) + FySf2.*repmat(nyJ(:,e)',Nq,1);
    FSf3 = FxSf3.*repmat(nxJ(:,e)',Nq,1) + FySf3.*repmat(nyJ(:,e)',Nq,1);    
    
    % bulk of GPU work: application of local operators
    fSproj1 = sum(VfPq.*FSf1',2);
    fSproj2 = sum(VfPq.*FSf2',2);
    fSproj3 = sum(VfPq.*FSf3',2);     
    
    % intermediate fluxes - form on the fly on GPU
    f1 = FSf1 + repmat((fSf1(:,e) - fSproj1)',Nq,1);
    f2 = FSf2 + repmat((fSf2(:,e) - fSproj2)',Nq,1);
    f3 = FSf3 + repmat((fSf3(:,e) - fSproj3)',Nq,1);
    
    %% project back to polynomial space - can incorporate WADG here
    
    JK = J(:,e);
    rhs1(:,e) =  Pq*((divF1 + .5*sum(VqLq.*f1,2))./JK);
    rhs2(:,e) =  Pq*((divF2 + .5*sum(VqLq.*f2,2))./JK);
    rhs3(:,e) =  Pq*((divF3 + .5*sum(VqLq.*f3,2))./JK);

end

% Lax-Friedrichs flux
g = 1;
cvel = sqrt(g*hM);
LFc = max(abs(sqrt(uM.^2+vM.^2))+cvel);

tau = 1;
hUn = (hP.*uP-hM.*uM).*nx + (hP.*vP-hM.*vM).*ny;
rhs1 = -2*rhs1 + .5*tau*Lq*(LFc.*(hP-hM).*Fscale);
rhs2 = -2*rhs2 + .5*tau*Lq*(LFc.*(hUn.*nx).*Fscale);
rhs3 = -2*rhs3 + .5*tau*Lq*(LFc.*(hUn.*ny).*Fscale);

end

