clear
Globals2D

N = 4;
K1D = 8;
FinalTime = .2;
CFL = .75;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
StartUp2D;

BuildPeriodicMaps2D(2,2);
% StartUp2D;

% plotting nodes
[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

global Vq Pq Lq Lqf Vfqf Vfq Pfqf
global rxJ sxJ ryJ syJ nxJ nyJ nrJ nsJ
global mapPq
Nq = 2*N+2;
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

%% problem params setup

x0 = 0; y0 = 0;

hex = @(x,y) 2+exp(-5^2*((x-x0).^2 + (y-y0).^2));
hex = @(x,y) 2 + (abs(x)<.5).*(abs(y)<.5);
h = Pq*hex(xq,yq);
hu = Pq*(xq*0);
hv = Pq*(xq*0);

% bathymetry + gravity
b = 0;
g = 1;

global fxS1 fxS2 fxS3 fyS1 fyS2 fyS3
avg = @(x,y) .5*(x+y);
fxS1 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(uL,uR);
fxS2 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(uL,uR).^2 + .5*g*avg(hL.^2,hR.^2);
fxS3 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(uL,uR).*avg(vL,vR);

fyS1 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(vL,vR);
fyS2 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(uL,uR).*avg(vL,vR);
fyS3 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(vL,vR).^2 + .5*g*avg(hL.^2,hR.^2);

%%

global wJq
wJq = diag(wq)*(Vq*J);

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
        
        hq = Vq*h;
        huq = Vq*hu;
        hvq = Vq*hv;
        uq = huq./hq;
        vq = hvq./hq;
        
        % project to entropy variables
        q1 = Pq*(g.*(hq+b)- .5*(uq.^2 + vq.^2));
        u  = Pq*uq;
        v  = Pq*vq;
        
        % re-evaluate at quad/surface points
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
    
    if mod(i,10)==0 || i==Nsteps
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
global rxJ sxJ ryJ syJ nxJ nyJ nrJ nsJ
global mapPq
global fxS1 fxS2 fxS3 fyS1 fyS2 fyS3

% operators
% Pq
% Vq
% Vfq Lq
Drq = (Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vfq*Pq);
Dsq = (Vq*Ds*Pq - .5*Vq*Lq*diag(nsJ)*Vfq*Pq);
% Lrq = (.5*Vq*Lq*diag(nrJ));
% Lsq = (.5*Vq*Lq*diag(nsJ));
VfPq = (Vfq*Pq);
VqLq = Vq*Lq;

hP = hM(mapPq);
uP = uM(mapPq);
vP = vM(mapPq);

fSf1 = nxJ.*fxS1(hM,uM,vM,hP,uP,vP) + nyJ.*fyS1(hM,uM,vM,hP,uP,vP);
fSf2 = nxJ.*fxS2(hM,uM,vM,hP,uP,vP) + nyJ.*fyS2(hM,uM,vM,hP,uP,vP);
fSf3 = nxJ.*fxS3(hM,uM,vM,hP,uP,vP) + nyJ.*fyS3(hM,uM,vM,hP,uP,vP);

rhs1 = zeros(Np,K);
rhs2 = zeros(Np,K);
rhs3 = zeros(Np,K);
for e = 1:K
    [hx hy] = meshgrid(hq(:,e));  [hfx hfy] = meshgrid(hM(:,e),hq(:,e));
    [ux uy] = meshgrid(uq(:,e));  [ufx ufy] = meshgrid(uM(:,e),uq(:,e));
    [vx vy] = meshgrid(vq(:,e));  [vfx vfy] = meshgrid(vM(:,e),vq(:,e));
    
    FxS1 = fxS1(hx,ux,vx,hy,uy,vy);  FyS1 = fyS1(hx,ux,vx,hy,uy,vy);      
    FxS2 = fxS2(hx,ux,vx,hy,uy,vy);  FyS2 = fyS2(hx,ux,vx,hy,uy,vy);  
    FxS3 = fxS3(hx,ux,vx,hy,uy,vy);  FyS3 = fyS3(hx,ux,vx,hy,uy,vy);  
    
    FxSf1 = fxS1(hfx,ufx,vfx,hfy,ufy,vfy);  FySf1 = fyS1(hfx,ufx,vfx,hfy,ufy,vfy);
    FxSf2 = fxS2(hfx,ufx,vfx,hfy,ufy,vfy);  FySf2 = fyS2(hfx,ufx,vfx,hfy,ufy,vfy);
    FxSf3 = fxS3(hfx,ufx,vfx,hfy,ufy,vfy);  FySf3 = fyS3(hfx,ufx,vfx,hfy,ufy,vfy);
       
    Nq = size(FxSf1,1);
    FSf1 = FxSf1.*repmat(nxJ(:,e)',Nq,1) + FySf1.*repmat(nyJ(:,e)',Nq,1);
    FSf2 = FxSf2.*repmat(nxJ(:,e)',Nq,1) + FySf2.*repmat(nyJ(:,e)',Nq,1);
    FSf3 = FxSf3.*repmat(nxJ(:,e)',Nq,1) + FySf3.*repmat(nyJ(:,e)',Nq,1);
        
%     keyboard
    
%     fx1r = sum(Drq.*FxS1,2) + sum(Lrq.*(FxSf1),2);  fx1s = sum(Dsq.*FxS1,2) + sum(Lsq.*(FxSf1),2);
%     fx2r = sum(Drq.*FxS2,2) + sum(Lrq.*(FxSf2),2);  fx2s = sum(Dsq.*FxS2,2) + sum(Lsq.*(FxSf2),2);
%     fx3r = sum(Drq.*FxS3,2) + sum(Lrq.*(FxSf3),2);  fx3s = sum(Dsq.*FxS3,2) + sum(Lsq.*(FxSf3),2);
%     
%     fy1r = sum(Drq.*FyS1,2) + sum(Lrq.*(FySf1),2);  fy1s = sum(Dsq.*FyS1,2) + sum(Lsq.*(FySf1),2);
%     fy2r = sum(Drq.*FyS2,2) + sum(Lrq.*(FySf2),2);  fy2s = sum(Dsq.*FyS2,2) + sum(Lsq.*(FySf2),2);
%     fy3r = sum(Drq.*FyS3,2) + sum(Lrq.*(FySf3),2);  fy3s = sum(Dsq.*FyS3,2) + sum(Lsq.*(FySf3),2);
         
    fx1r = sum(Drq.*FxS1,2);  fx1s = sum(Dsq.*FxS1,2);
    fx2r = sum(Drq.*FxS2,2);  fx2s = sum(Dsq.*FxS2,2);
    fx3r = sum(Drq.*FxS3,2);  fx3s = sum(Dsq.*FxS3,2);
    
    fy1r = sum(Drq.*FyS1,2);  fy1s = sum(Dsq.*FyS1,2);
    fy2r = sum(Drq.*FyS2,2);  fy2s = sum(Dsq.*FyS2,2);
    fy3r = sum(Drq.*FyS3,2);  fy3s = sum(Dsq.*FyS3,2);

    % apply J*G^T here
    divF1 = rxJ(:,e).*fx1r + sxJ(:,e).*fx1s + ryJ(:,e).*fy1r + syJ(:,e).*fy1s;
    divF2 = rxJ(:,e).*fx2r + sxJ(:,e).*fx2s + ryJ(:,e).*fy2r + syJ(:,e).*fy2s;
    divF3 = rxJ(:,e).*fx3r + sxJ(:,e).*fx3s + ryJ(:,e).*fy3r + syJ(:,e).*fy3s;
    
    %fSproj1 = nxJ(:,e).*sum(VfPq.*FxSf1',2) + nyJ(:,e).*sum(VfPq.*FySf1',2);
    %fSproj2 = nxJ(:,e).*sum(VfPq.*FxSf2',2) + nyJ(:,e).*sum(VfPq.*FySf2',2);
    %fSproj3 = nxJ(:,e).*sum(VfPq.*FxSf3',2) + nyJ(:,e).*sum(VfPq.*FySf3',2);
    fSproj1 = sum(VfPq.*FSf1',2);
    fSproj2 = sum(VfPq.*FSf2',2);
    fSproj3 = sum(VfPq.*FSf3',2);
    
%     rhs1(:,e) =  Pq*divF1 + Lq*(.5*(fSf1(:,e) - fSproj1));
%     rhs2(:,e) =  Pq*divF2 + Lq*(.5*(fSf2(:,e) - fSproj2));
%     rhs3(:,e) =  Pq*divF3 + Lq*(.5*(fSf3(:,e) - fSproj3));
        
    f1 = .5*FSf1 + repmat(.5*(fSf1(:,e) - fSproj1)',Nq,1);
    f2 = .5*FSf2 + repmat(.5*(fSf2(:,e) - fSproj2)',Nq,1);
    f3 = .5*FSf3 + repmat(.5*(fSf3(:,e) - fSproj3)',Nq,1);
    
    rhs1(:,e) =  Pq*(divF1 + sum(VqLq.*f1,2));
    rhs2(:,e) =  Pq*(divF2 + sum(VqLq.*f2,2));
    rhs3(:,e) =  Pq*(divF3 + sum(VqLq.*f3,2));

end

% Lax-Friedrichs flux
g = 1;
cvel = sqrt(g*hM);
LFc = max(abs(uM)+cvel);

tau = 1;
hUn = (hP.*uP-hM.*uM).*nx + (hP.*vP-hM.*vM).*ny;
rhs1 = -2*rhs1./J + .5*tau*Lq*(LFc.*(hP-hM).*Fscale);
rhs2 = -2*rhs2./J + .5*tau*Lq*(LFc.*(hUn.*nx).*Fscale);
rhs3 = -2*rhs3./J + .5*tau*Lq*(LFc.*(hUn.*ny).*Fscale);

end

