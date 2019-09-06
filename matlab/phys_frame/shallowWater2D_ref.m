clear
Globals2D

a = 0/8; % warping factor
% FinalTime = 10;

N = 3;
K1D = 8;
FinalTime = 1;
plotMesh = 0;

wadgProjEntropyVars = abs(a)>1e-8;
CFL = .5;
global tau
tau = 1;

Lx = 1; Ly = 1; ratiox = 1; ratioy = 1;
% Lx = 7.5; Ly = 5; ratiox = 3/4; ratioy = .5;
% Lx = 10; Ly = 5; ratiox = 1; ratioy = Ly/Lx;
[Nv, VX, VY, K, EToV] = unif_tri_mesh(round(ratiox*K1D),round(K1D*ratioy));
VX = VX/max(abs(VX));  VY = VY/max(abs(VY));
VX = (VX+1)*Lx; VY = VY*Ly;

% VX = VX + randn(size(VX));
StartUp2D;
BuildPeriodicMaps2D(max(VX)-min(VX),max(VY)-min(VY));

%% reference operators

% quadrature nodes
Nq = 2*N;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
[rq1D wq1D] = JacobiGQ(0,0,N);
Vq = Vandermonde2D(N,rq,sq)/V;

% plotting nodes
[rp sp] = EquiNodes2D(15); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;

rfq = [rq1D; -rq1D; -ones(size(rq1D))];
sfq = [-ones(size(rq1D)); rq1D; -rq1D];
wfq = [wq1D; wq1D; wq1D];
Vq1D = Vandermonde1D(N,rq1D)/Vandermonde1D(N,JacobiGL(0,0,N));

Vff = kron(eye(3),Vq1D);
Vf = Vandermonde2D(N,rfq,sfq)/V;

nrJ = [-zeros(size(rq1D)); ones(size(rq1D)); -ones(size(rq1D)); ];
nsJ = [-ones(size(rq1D)); ones(size(rq1D)); -zeros(size(rq1D)); ];

Nfq = length(rq1D);
Nq = length(rq);

%% flux differencing operators

M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq)); % J's cancel out
Lq = M\(Vf'*diag(wfq));

Qr = diag(wq)*Vq*Dr*Pq;
Qs = diag(wq)*Vq*Ds*Pq;
Br = diag(nrJ.*wfq);
Bs = diag(nsJ.*wfq);
E = (Vf*Pq);

global QNr QNs VqPN VqPq VqLq
QNr = [Qr-.5*E'*Br*E .5*E'*Br;
    -.5*Br*E .5*Br];
QNs = [Qs-.5*E'*Bs*E .5*E'*Bs;
    -.5*Bs*E .5*Bs];

% make skew symmetric operators
DNr = diag(1./[wq;wfq])*(QNr-QNr')*.5;
DNs = diag(1./[wq;wfq])*(QNs-QNs')*.5;

% projection operators
VqPq = Vq*Pq;
VqLq = Vq*Lq;
VqPN = [VqPq VqLq];

%% make quadrature face maps

global mapPq

% get face nodes
xf = Vf*x;
yf = Vf*y;

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

% optional - find wall boundary nodes
global mapBwall
%mapBwall = find(abs(xf-max(x(:)))<1e-8 | abs(xf-min(x(:)))<1e-8);
% mapBwall = find(abs(yf-max(y(:)))<1e-8 | abs(yf-min(y(:)))<1e-8);
mapBwall = [];

% plot(xf,yf,'o')
% hold on
% plot(xf(mapBwall),yf(mapBwall),'x','markersize',15)
% return


%% make curvilinear mesh 

x0 = Lx; y0 = 0;
x = x + Lx*a*cos(1/2*pi*(x-x0)/Lx).*cos(3/2*pi*(y-y0)/Ly);
y = y + Ly*a*sin(pi*(x-x0)/Lx).*cos(1/2*pi*(y-y0)/Ly);

xq = Vq*x; yq = Vq*y;
xp = Vp*x; yp = Vp*y;
xf = Vf*x; yf = Vf*y;

if plotMesh
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

% curved geometric terms
global rxJ sxJ ryJ syJ nxJ nyJ 
[rxk,sxk,ryk,syk,Jk] = GeometricFactors2D(x,y,[Vq;Vf]*Dr,[Vq;Vf]*Ds);
rxJ = rxk.*Jk;    sxJ = sxk.*Jk;
ryJ = ryk.*Jk;    syJ = syk.*Jk;
J = Jk(1:Nq,:);
Jf = Jk(Nq+1:end,:);
nxJ = rxJ(Nq+1:end,:).*nrJ + sxJ(Nq+1:end,:).*nsJ;
nyJ = ryJ(Nq+1:end,:).*nrJ + syJ(Nq+1:end,:).*nsJ;
sJ = sqrt(nxJ.^2 + nyJ.^2);
nx = nxJ./sJ; ny = nyJ./sJ;


%% entropy stable fluxes

global g b
g = 1;
b = 0;

% entropy variables
V1 = @(h,hu,hv) g.*(h+b)- .5*((hu./h).^2 + (hv./h).^2);
V2 = @(h,hu,hv) hu./h;
V3 = @(h,hu,hv) hv./h;
U1 = @(v1,v2,v3) (v1+.5*(v2.^2+v3.^2))/g - b;
U2 = @(v1,v2,v3) U1(v1,v2,v3).*v2;
U3 = @(v1,v2,v3) U1(v1,v2,v3).*v3;

global fxS1 fyS1 fxS2 fyS2 fxS3 fyS3 avg
avg = @(x,y) .5*(x+y);
fxS1 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(uL,uR);
fxS2 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(uL,uR).^2 + .5*g*avg(hL.^2,hR.^2);
fxS3 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(uL,uR).*avg(vL,vR);

fyS1 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(vL,vR);
fyS2 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(uL,uR).*avg(vL,vR);
fyS3 = @(hL,uL,vL,hR,uR,vR) avg(hL,hR).*avg(vL,vR).^2 + .5*g*avg(hL.^2,hR.^2);

% initial condition
h = 2+exp(-50*((xq-x0).^2+(yq-y0).^2));
hu = 0*h;
hv = 0*h;

%%

global wJq
wJq = diag(wq)*(J);

% Runge-Kutta residual storage
res1 = zeros(Nq,K);
res2 = zeros(Nq,K);
res3 = zeros(Nq,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
CNh = max(CN*max(sJ(:)./Jf(:)));
dt = CFL*2/CNh;

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

figure(1)
for i = 1:Nsteps
    for INTRK = 1:5               
        
        % project to entropy variables for curvilinear
        if wadgProjEntropyVars
            q1 = Pq*((VqPq*(V1(h,hu,hv).*J))./J);
            q2 = Pq*((VqPq*(V2(h,hu,hv).*J))./J);
            q3 = Pq*((VqPq*(V3(h,hu,hv).*J))./J);
        else
            q1 = Pq*V1(h,hu,hv);
            q2 = Pq*V2(h,hu,hv);
            q3 = Pq*V3(h,hu,hv);            
        end
        
        % evaluate at quad/surface points
        q1q = Vq*q1;  q1M = Vf*q1;
        q2q = Vq*q2;  q2M = Vf*q2;
        q3q = Vq*q3;  q3M = Vf*q3;        
        hq  = U1(q1q,q2q,q3q); 
        huq = U2(q1q,q2q,q3q); 
        hvq = U3(q1q,q2q,q3q); 
        hM  = U1(q1M,q2M,q3M);
        huM = U2(q1M,q2M,q3M);
        hvM = U3(q1M,q2M,q3M);        
        
        uq = huq./hq; uM = huM./hM;
        vq = hvq./hq; vM = hvM./hM;
        
        % extra LF flux info        
        [rhs1 rhs2 rhs3]  = RHS2Dsimple(hq,uq,vq,hM,uM,vM);
        
        if (INTRK==5)
            rhstest(i) = 0;
            for e = 1:K                
                q1e = J(:,e).*wq.*(Vq*q1(:,e));
                q2e = J(:,e).*wq.*(Vq*q2(:,e));
                q3e = J(:,e).*wq.*(Vq*q3(:,e));                
                r1 = rhs1(:,e);
                r2 = rhs2(:,e);
                r3 = rhs3(:,e);                

                rhstest(i) = rhstest(i) + (q1e'*r1 + q2e'*r2 + q3e'*r3);
            end
        end
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;        
        
        h  = h  + rk4b(INTRK)*res1;
        hu = hu + rk4b(INTRK)*res2;
        hv = hv + rk4b(INTRK)*res3;        
        
    end;
        
    if  (mod(i,5)==0 || i==Nsteps)
        clf
        pp = h;
        vv = real(Vp*Pq*pp);
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f, N = %d, K1D = %d',dt*i,N,K1D))
%         view(3)
        drawnow
    end    
end

figure(2)
semilogy(dt*(1:Nsteps),abs(rhstest),'x--');hold on
legend('rhstest')


function [rhs1 rhs2 rhs3 rhs4] = RHS2Dsimple(hq,uq,vq,hM,uM,vM)

Globals2D;

global rxJ sxJ ryJ syJ nxJ nyJ 
global mapPq
global fxS1 fyS1 fxS2 fyS2 fxS3 fyS3 avg
global DNr DNs VqPq VqLq VqPN

hP = hM(mapPq);
uP = uM(mapPq);
vP = vM(mapPq);

global g b tau
LFc = sqrt(g*hM) + abs(uM.*nx + vM.*ny); % max wavespeed
dQ1 = hP-hM;
dQ2 = hP.*uP-hM.*uM;
dQ3 = hP.*vP-hM.*vM;
Lf1 = tau*LFc.*dQ1.*sJ;
Lf2 = tau*LFc.*dQ2.*sJ;
Lf3 = tau*LFc.*dQ3.*sJ;

fSf1 = nxJ.*fxS1(hM,uM,vM,hP,uP,vP) + nyJ.*fyS1(hM,uM,vM,hP,uP,vP);
fSf2 = nxJ.*fxS2(hM,uM,vM,hP,uP,vP) + nyJ.*fyS2(hM,uM,vM,hP,uP,vP);
fSf3 = nxJ.*fxS3(hM,uM,vM,hP,uP,vP) + nyJ.*fyS3(hM,uM,vM,hP,uP,vP);

fSf1 = fSf1  - .25*Lf1;
fSf2 = fSf2  - .25*Lf2;
fSf3 = fSf3  - .25*Lf3;

h = [hq; hM];
u = [uq; uM];
v = [vq; vM];
divF1 = zeros(size(DNr,1),K);
divF2 = zeros(size(DNr,1),K);
divF3 = zeros(size(DNr,1),K);
for e = 1:K
    
    [hx, hy] = meshgrid(h(:,e));
    [ux, uy] = meshgrid(u(:,e));
    [vx, vy] = meshgrid(v(:,e));
    
    % optimized evaluations
    ha = avg(hx,hy);
    h2a = avg(hx.^2,hy.^2);
    ua = avg(ux,uy);
    va = avg(vx,vy);
    FxS1 = ha.*ua;
    FxS2 = FxS1.*ua + .5*g*h2a;
    FxS3 = FxS1.*va;    
    
    FyS1 = ha.*va;
    FyS2 = FyS1.*ua;
    FyS3 = FyS1.*va + .5*g*h2a;
        
    [rxJ1, rxJ2] = meshgrid(rxJ(:,e));
    [sxJ1, sxJ2] = meshgrid(sxJ(:,e));
    [ryJ1, ryJ2] = meshgrid(ryJ(:,e));
    [syJ1, syJ2] = meshgrid(syJ(:,e));
    rxJK = avg(rxJ1,rxJ2);  sxJK = avg(sxJ1,sxJ2);
    ryJK = avg(ryJ1,ryJ2);  syJK = avg(syJ1,syJ2);
    
    Dx = DNr.*rxJK + DNs.*sxJK;
    Dy = DNr.*ryJK + DNs.*syJK;
    
    divF1(:,e) = sum(Dx.*FxS1,2) + sum(Dy.*FyS1,2);
    divF2(:,e) = sum(Dx.*FxS2,2) + sum(Dy.*FyS2,2);
    divF3(:,e) = sum(Dx.*FxS3,2) + sum(Dy.*FyS3,2);   

end

rhs1 = 2*VqPN*divF1 + VqLq*(fSf1);
rhs2 = 2*VqPN*divF2 + VqLq*(fSf2);
rhs3 = 2*VqPN*divF3 + VqLq*(fSf3);

% apply wadg
rhs1 = -VqPq*(rhs1./J);
rhs2 = -VqPq*(rhs2./J);
rhs3 = -VqPq*(rhs3./J);


end
