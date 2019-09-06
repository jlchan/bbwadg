% check_esdg_3d
clear
clear -globals

Globals3D

N = 1;
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGmsh3D('~/Desktop/wadges/meshes/periodicCubeTGCoarse.msh');
StartUp3D;

[rq sq tq wq] = tet_cubature(2*N);
Nq = length(rq);
Vq = Vandermonde3D(N,rq,sq,tq)/V;
M =(Vq'*diag(wq)*Vq);
Pq = M\(Vq'*diag(wq));

xq = Vq*x;
yq = Vq*y;
zq = Vq*z;

Drq = Vq*Dr;
Dsq = Vq*Ds;
Dtq = Vq*Dt;

Nfp = (N+1)*(N+2)/2;
[rqtri sqtri wqtri] = Cubature2D(2*N);
Vfqf = Vandermonde2D(N,rqtri,sqtri)/Vandermonde2D(N,r(Fmask(:,1)),s(Fmask(:,1)));
Vfqf = kron(eye(4),Vfqf);
rfq = Vfqf*r(Fmask(:));
sfq = Vfqf*s(Fmask(:));
tfq = Vfqf*t(Fmask(:));
wfq = repmat(wqtri,4,1);
Vfq = Vandermonde3D(N,rfq,sfq,tfq)/V;
Lq = M\(Vfq'*diag(wfq));
VqPq = Vq*Pq;
VqLq = Vq*Lq;
VfPq = Vfq*Pq;

Drfq = Vfq*Dr;
Dsfq = Vfq*Ds;
Dtfq = Vfq*Dt;

Nfqf = length(rqtri);
e = ones(Nfqf,1); zz = zeros(Nfqf,1);
nrJ = [zz;zz;e;-e];
nsJ = [zz;-e;e;zz];
ntJ = [-e;zz;e;zz];

DNr = [Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vfq*Pq .5*Vq*Lq*diag(nrJ);
    -.5*diag(nrJ)*Vfq*Pq .5*diag(nrJ)];
DNs = [Vq*Ds*Pq - .5*Vq*Lq*diag(nsJ)*Vfq*Pq .5*Vq*Lq*diag(nsJ);
    -.5*diag(nsJ)*Vfq*Pq .5*diag(nsJ)];
DNt = [Vq*Dt*Pq - .5*Vq*Lq*diag(ntJ)*Vfq*Pq .5*Vq*Lq*diag(ntJ);
    -.5*diag(ntJ)*Vfq*Pq .5*diag(ntJ)];
WN = diag([wq;wfq]);

DNr_skew = [Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vfq*Pq .5*Vq*Lq*diag(nrJ);
    -.5*diag(nrJ)*Vfq*Pq 0*.5*diag(nrJ)];
QNr_skew = WN*DNr_skew;

QNr = WN*DNr;
QNs = WN*DNs;
QNt = WN*DNt;

% DNrq = Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vfq*Pq
% DNsq = Vq*Ds*Pq - .5*Vq*Lq*diag(nsJ)*Vfq*Pq
% DNtq = Vq*Dt*Pq - .5*Vq*Lq*diag(ntJ)*Vfq*Pq

rxJf = Vfq*(rx.*J); sxJf = Vfq*(sx.*J); txJf = Vfq*(tx.*J);
ryJf = Vfq*(ry.*J); syJf = Vfq*(sy.*J); tyJf = Vfq*(ty.*J);
rzJf = Vfq*(rz.*J); szJf = Vfq*(sz.*J); tzJf = Vfq*(tz.*J);
nxJ = rxJf.*nrJ + sxJf.*nsJ + txJf.*ntJ;
nyJ = ryJf.*nrJ + syJf.*nsJ + tyJf.*ntJ;
nzJ = rzJf.*nrJ + szJf.*nsJ + tzJf.*ntJ;

rxJ = Vq*(rx.*J); sxJ = Vq*(sx.*J); txJ = Vq*(tx.*J);
ryJ = Vq*(ry.*J); syJ = Vq*(sy.*J); tyJ = Vq*(ty.*J);
rzJ = Vq*(rz.*J); szJ = Vq*(sz.*J); tzJ = Vq*(tz.*J);

DNrmod = [Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vfq*Pq .5*Vq*Lq*diag(nrJ);
    -.5*diag(nrJ)*Vfq*Pq 0*.5*diag(nrJ)];
DNsmod = [Vq*Ds*Pq - .5*Vq*Lq*diag(nsJ)*Vfq*Pq .5*Vq*Lq*diag(nsJ);
    -.5*diag(nsJ)*Vfq*Pq 0*.5*diag(nsJ)];
DNtmod = [Vq*Dt*Pq - .5*Vq*Lq*diag(ntJ)*Vfq*Pq .5*Vq*Lq*diag(ntJ);
    -.5*diag(ntJ)*Vfq*Pq 0*.5*diag(ntJ)];

%% make quadrature face maps

Nfq = size(Vfq,1)/Nfaces;
xf = Vfq*x;
yf = Vfq*y;
zf = Vfq*z;

mapMq = reshape(1:length(xf(:)),Nfq*Nfaces,K);
mapPq = mapMq;

for e = 1:K
    for f = 1:Nfaces
        enbr = EToE(e,f);
        if e ~= enbr % if it's a neighbor
            fnbr = EToF(e,f);
            id1 = (1:Nfq) + (f-1)*Nfq;
            id2 = (1:Nfq) + (fnbr-1)*Nfq;
            x1 = xf(id1,e); y1 = yf(id1,e); z1 = zf(id1,e);
            x2 = xf(id2,enbr); y2 = yf(id2,enbr); z2 = zf(id2,enbr);
            
            % find matches
            p = zeros(Nfq,1);
            for i = 1:Nfq
                for j = 1:Nfq
                    d = abs(x1(i)-x2(j)) + abs(y1(i)-y2(j)) + abs(z1(i)-z2(j));
                    if d < 1e-8
                        p(i) = j;
                    end
                end
            end
            mapPq(id1,e) = id2(p) + (enbr-1)*(Nfq*Nfaces);
        else
            % if e is on boundary
            
        end
        
        
    end
end
mapBq = find(mapPq==mapMq);

% norm(xf(mapMq)-xf(mapPq),'fro')
% norm(yf(mapMq)-yf(mapPq),'fro')
% norm(zf(mapMq)-zf(mapPq),'fro')

%%

clear mapM mapP
esdg_data;
mapMq = mapM;
mapPq = mapP;

%%

% d = (xf-xf(mapPq)).^2 + (yf-yf(mapPq)).^2 + (zf-zf(mapPq)).^2;
% ids = find(d > 1e-12);
% mapPq(ids) = mapMq(ids);

% bids = ids;
% bids = find(((abs(xf)<1e-8)) | (abs(xf-10)<1e-8) | abs(yf)<1e-8 | abs(yf-20)<1e-8 | abs(zf)<1e-8 | abs(zf-10)<1e-8);
% mapPq(bids)=mapMq(bids);
% bids = find(mapPq==mapMq);
% plot3(xf(bids),yf(bids),zf(bids),'o')

d = 0;
for i = 1:Nfq*Nfaces*K
    idM = mapMq(i);
    idP = mapPq(i);
    if (abs(nxJ(idM)+nxJ(idP)) + abs(nyJ(idM)+nyJ(idP))+ abs(nzJ(idM)+nzJ(idP))) > 1e-11 && idM~=idP
        [nxJ(idM),nxJ(idP)]
        [nyJ(idM),nyJ(idP)]
        [nzJ(idM),nzJ(idP)]
        keyboard
    end
    
%     clf
%     plot3(xf(idM),yf(idM),zf(idM),'o')
%     hold on
%     plot3(xf(idP),yf(idP),zf(idP),'x')
%     quiver3(xf(idM),yf(idM),zf(idM),nxJ(idM),nyJ(idM),nzJ(idM))
%     quiver3(xf(idP),yf(idP),zf(idP),nxJ(idP),nyJ(idP),nzJ(idP))
%     title(['idM = ' num2str(i)] )
%     axis([-1 1 -1 1 -1 1]*pi*1.1)
% %     axis equal
%     pause(.1)
end
% return


%% check 3D entropy fluxes

global gamma avg pfun betafun 
gamma = 1.4;
avg = @(uL,uR) .5*(uL+uR);
rhoe = @(rho,rhou,rhov,rhow,E) E - .5*(rhou.^2+rhov.^2+rhow.^2)./rho;
pcons = @(rho,rhou,rhov,rhow,E) (gamma-1)*rhoe(rho,rhou,rhov,rhow,E);
s = @(rho,rhou,rhov,rhow,E) log((gamma-1)*rhoe(rho,rhou,rhov,rhow,E)./(rho.^gamma));
pfun = @(rho, u, v, w, E) ((gamma - 1.0) * (E - .5 * rho .* (u .* u + v .* v + w .* w)));
betafun = @(rho, u, v, w, E) (rho ./ (2.0 * pfun(rho, u, v, w, E)));

sV = @(V1,V2,V3,V4,V5) gamma - V1 + (V2.^2+V3.^2+V4.^2)./(2*V5);
rhoeV  = @(V1,V2,V3,V4,V5) ((gamma-1)./((-V5).^gamma)).^(1/(gamma-1)).*exp(-sV(V1,V2,V3,V4,V5)/(gamma-1));

V1 = @(rho,rhou,rhov,rhow,E) (-E + rhoe(rho,rhou,rhov,rhow,E).*(gamma + 1 - s(rho,rhou,rhov,rhow,E)))./(rhoe(rho,rhou,rhov,rhow,E));
V2 = @(rho,rhou,rhov,rhow,E) rhou./(rhoe(rho,rhou,rhov,rhow,E));
V3 = @(rho,rhou,rhov,rhow,E) rhov./(rhoe(rho,rhou,rhov,rhow,E));
V4 = @(rho,rhou,rhov,rhow,E) rhow./(rhoe(rho,rhou,rhov,rhow,E));
V5 = @(rho,rhou,rhov,rhow,E) (-rho)./(rhoe(rho,rhou,rhov,rhow,E));

U1 = @(V1,V2,V3,V4,V5) rhoeV(V1,V2,V3,V4,V5).*(-V5);
U2 = @(V1,V2,V3,V4,V5) rhoeV(V1,V2,V3,V4,V5).*(V2);
U3 = @(V1,V2,V3,V4,V5) rhoeV(V1,V2,V3,V4,V5).*(V3);
U4 = @(V1,V2,V3,V4,V5) rhoeV(V1,V2,V3,V4,V5).*(V4);
U5 = @(V1,V2,V3,V4,V5) rhoeV(V1,V2,V3,V4,V5).*(1-(V2.^2+V3.^2+V4.^2)./(2*V5));

% entropy potentials
psix = @(rho,rhou,rhov,rhow,E) (gamma-1)*rhou;
psiy = @(rho,rhou,rhov,rhow,E) (gamma-1)*rhov;
psiz = @(rho,rhou,rhov,rhow,E) (gamma-1)*rhow;

%% test 2 point flux condition

if 1
    rhoL = 2+rand;
    uL = randn;
    vL = randn;
    wL = randn;
    EL = 1 + .5*rhoL.*(uL.^2+vL.^2+wL.^2);
    rhoR = 2;
    uR = 3;
    vR = 4;
    wR = 5;
    ER = 1 + .5*rhoR.*(uR.^2+vR.^2+wR.^2);
    
    rhouL = rhoL*uL;
    rhovL = rhoL*vL;
    rhowL = rhoL*wL;
    rhouR = rhoR*uR;
    rhovR = rhoR*vR;
    rhowR = rhoR*wR;
    
    UVU = [U1(V1(rhoL,rhouL,rhovL,rhowL,EL),V2(rhoL,rhouL,rhovL,rhowL,EL),V3(rhoL,rhouL,rhovL,rhowL,EL),V4(rhoL,rhouL,rhovL,rhowL,EL),V5(rhoL,rhouL,rhovL,rhowL,EL))
        U2(V1(rhoL,rhouL,rhovL,rhowL,EL),V2(rhoL,rhouL,rhovL,rhowL,EL),V3(rhoL,rhouL,rhovL,rhowL,EL),V4(rhoL,rhouL,rhovL,rhowL,EL),V5(rhoL,rhouL,rhovL,rhowL,EL))
        U3(V1(rhoL,rhouL,rhovL,rhowL,EL),V2(rhoL,rhouL,rhovL,rhowL,EL),V3(rhoL,rhouL,rhovL,rhowL,EL),V4(rhoL,rhouL,rhovL,rhowL,EL),V5(rhoL,rhouL,rhovL,rhowL,EL))
        U4(V1(rhoL,rhouL,rhovL,rhowL,EL),V2(rhoL,rhouL,rhovL,rhowL,EL),V3(rhoL,rhouL,rhovL,rhowL,EL),V4(rhoL,rhouL,rhovL,rhowL,EL),V5(rhoL,rhouL,rhovL,rhowL,EL))
        U5(V1(rhoL,rhouL,rhovL,rhowL,EL),V2(rhoL,rhouL,rhovL,rhowL,EL),V3(rhoL,rhouL,rhovL,rhowL,EL),V4(rhoL,rhouL,rhovL,rhowL,EL),V5(rhoL,rhouL,rhovL,rhowL,EL))];
    
    UU = [rhoL;rhouL;rhovL;rhowL;EL];
    
    norm(UVU-UU)
    
    [FxS1 FxS2 FxS3 FxS4 FxS5 FyS1 FyS2 FyS3 FyS4 FyS5 FzS1 FzS2 FzS3 FzS4 FzS5] = fluxes(rhoL,rhoR,uL,uR,vL,vR,wL,wR,EL,ER);
    
    V1L = V1(rhoL,rhoL*uL,rhoL*vL,rhoL*wL,EL);  V1R = V1(rhoR,rhoR*uR,rhoR*vR,rhoR*wR,ER);
    V2L = V2(rhoL,rhoL*uL,rhoL*vL,rhoL*wL,EL);  V2R = V2(rhoR,rhoR*uR,rhoR*vR,rhoR*wR,ER);
    V3L = V3(rhoL,rhoL*uL,rhoL*vL,rhoL*wL,EL);  V3R = V3(rhoR,rhoR*uR,rhoR*vR,rhoR*wR,ER);
    V4L = V4(rhoL,rhoL*uL,rhoL*vL,rhoL*wL,EL);  V4R = V4(rhoR,rhoR*uR,rhoR*vR,rhoR*wR,ER);
    V5L = V5(rhoL,rhoL*uL,rhoL*vL,rhoL*wL,EL);  V5R = V5(rhoR,rhoR*uR,rhoR*vR,rhoR*wR,ER);
    
    a = (V1L-V1R).*FxS1 + (V2L-V2R).*FxS2 + (V3L-V3R).*FxS3 + (V4L-V4R).*FxS4 + (V5L-V5R).*FxS5;
    b = psix(rhoL,rhouL,rhovL,rhowL,EL)-psix(rhoR,rhouR,rhovR,rhowR,ER);
    abs(a-b)
    
    a = (V1L-V1R).*FyS1 + (V2L-V2R).*FyS2 + (V3L-V3R).*FyS3 + (V4L-V4R).*FyS4 + (V5L-V5R).*FyS5;
    b = psiy(rhoL,rhouL,rhovL,rhowL,EL)-psiy(rhoR,rhouR,rhovR,rhowR,ER);
    abs(a-b)
    
    a = (V1L-V1R).*FzS1 + (V2L-V2R).*FzS2 + (V3L-V3R).*FzS3 + (V4L-V4R).*FzS4 + (V5L-V5R).*FzS5;
    b = psiz(rhoL,rhouL,rhovL,rhowL,EL)-psiz(rhoR,rhouR,rhovR,rhowR,ER);
    abs(a-b)
    return
end


%% test RHS eval

rho = ones(size(xq));
u = sin(xq).*cos(yq).*cos(zq);
v = -cos(xq).*sin(yq).*cos(zq);
w = 0.0;
p = 100.0/gamma + (1.0/16.0).*(cos(2.0*xq) + cos(2.0*yq)).*(2.0 + cos(2.0*zq));

rho = 2+.0*rand(size(xq));
u   = 1+.0*rand(size(xq));
v   = 1+.0*rand(size(xq));
w   = 1+.1*rand(size(xq));
w(xq>0) = 1;
w(xq<0) = 2;
p   = 2+.0*rand(size(xq));

% x0 = pi/2; y0 = pi/2; z0 = pi/2;
% du = exp(-5^2*((xq-x0).^2+(yq-y0).^2+(zq-z0).^2));
% rho = 4 + du;
% u = 0 + du;
% v = 1 + du;
% w = 2 + du;
% p = 1 + du;

rho = Vq*Pq*rho;
rhou = Vq*Pq*(rho.*u);
rhov = Vq*Pq*(rho.*v);
rhow = Vq*Pq*(rho.*w);
E = Vq*Pq*(p/(gamma-1)+.5*rho.*(u.^2+v.^2+w.^2));

% rhof = Vfq*Pq*rho;
% rhouf = Vfq*Pq*rhou;
% rhovf = Vfq*Pq*rhov;
% rhowf = Vfq*Pq*rhow;
% Ef = Vfq*Pq*E;

%% test RHS eval

rho_o = rho;

v1 = VqPq*V1(rho,rhou,rhov,rhow,E);
v2 = VqPq*V2(rho,rhou,rhov,rhow,E);
v3 = VqPq*V3(rho,rhou,rhov,rhow,E);
v4 = VqPq*V4(rho,rhou,rhov,rhow,E);
v5 = VqPq*V5(rho,rhou,rhov,rhow,E);

rho = U1(v1,v2,v3,v4,v5);
rhou = U2(v1,v2,v3,v4,v5);
rhov = U3(v1,v2,v3,v4,v5);
rhow = U4(v1,v2,v3,v4,v5);
E = U5(v1,v2,v3,v4,v5);

u = rhou./rho;
v = rhov./rho;
w = rhow./rho;

rhof = U1(VfPq*v1,VfPq*v2,VfPq*v3,VfPq*v4,VfPq*v5);
rhouf = U2(VfPq*v1,VfPq*v2,VfPq*v3,VfPq*v4,VfPq*v5);
rhovf = U3(VfPq*v1,VfPq*v2,VfPq*v3,VfPq*v4,VfPq*v5);
rhowf = U4(VfPq*v1,VfPq*v2,VfPq*v3,VfPq*v4,VfPq*v5);
Ef = U5(VfPq*v1,VfPq*v2,VfPq*v3,VfPq*v4,VfPq*v5);

uf = rhouf./rhof;
vf = rhovf./rhof;
wf = rhowf./rhof;

% compute fluxes
rhoP = rhof(mapPq); rhoM = rhof(mapMq);
uP = uf(mapPq); uM = uf(mapMq);
vP = vf(mapPq); vM = vf(mapMq);
wP = wf(mapPq); wM = wf(mapMq);
EP = Ef(mapPq); EM = Ef(mapMq);

[FxSf1 FxSf2 FxSf3 FxSf4 FxSf5 FySf1 FySf2 FySf3 FySf4 FySf5 FzSf1 FzSf2 FzSf3 FzSf4 FzSf5] = fluxes(rhoM,rhoP,uM,uP,vM,vP,wM,wP,EM,EP);

r1 = zeros(size(DNr,1),K);
r2 = zeros(size(DNr,1),K);
r3 = zeros(size(DNr,1),K);
r4 = zeros(size(DNr,1),K);
r5 = zeros(size(DNr,1),K);
for e = 1:K
    
    [rhoL rhoR] = meshgrid([rho(:,e);rhof(:,e)]);
    [uL uR] = meshgrid([u(:,e);uf(:,e)]);
    [vL vR] = meshgrid([v(:,e);vf(:,e)]);
    [wL wR] = meshgrid([w(:,e);wf(:,e)]);
    [EL ER] = meshgrid([E(:,e);Ef(:,e)]);
    
    [FxS1 FxS2 FxS3 FxS4 FxS5 ...
        FyS1 FyS2 FyS3 FyS4 FyS5...
        FzS1 FzS2 FzS3 FzS4 FzS5] = fluxes(rhoL,rhoR,uL,uR,vL,vR,wL,wR,EL,ER);
    
    DNx = diag([rxJ(:,e);rxJf(:,e)])*DNrmod + diag([sxJ(:,e);sxJf(:,e)])*DNsmod + diag([txJ(:,e);txJf(:,e)])*DNtmod;
    DNy = diag([ryJ(:,e);ryJf(:,e)])*DNrmod + diag([syJ(:,e);syJf(:,e)])*DNsmod + diag([tyJ(:,e);tyJf(:,e)])*DNtmod;
    DNz = diag([rzJ(:,e);rzJf(:,e)])*DNrmod + diag([szJ(:,e);szJf(:,e)])*DNsmod + diag([tzJ(:,e);tzJf(:,e)])*DNtmod;
    
    r1(:,e) = sum(DNx.*FxS1,2) + sum(DNy.*FyS1,2) + sum(DNz.*FzS1,2);    
    r2(:,e) = sum(DNx.*FxS2,2) + sum(DNy.*FyS2,2) + sum(DNz.*FzS2,2);
    r3(:,e) = sum(DNx.*FxS3,2) + sum(DNy.*FyS3,2) + sum(DNz.*FzS3,2);
    r4(:,e) = sum(DNx.*FxS4,2) + sum(DNy.*FyS4,2) + sum(DNz.*FzS4,2);
    r5(:,e) = sum(DNx.*FxS5,2) + sum(DNy.*FyS5,2) + sum(DNz.*FzS5,2);    
    
%     FxS1 = ones(size(FxS1));
%     r1(:,e) = sum(DNr.*FxS1,2);
end

Nq = size(Vq,1);
% [sum(sum(r1(1:Nq,:))) sum(sum(r2(1:Nq,:))) sum(sum(r3(1:Nq,:))) sum(sum(r4(1:Nq,:))) sum(sum(r5(1:Nq,:)))]
% sum(sum(r1(Nq+1:end,:))) % face pts rhs

% should scale total by 2x
rhs1 = 2*([VqPq VqLq]*r1 + .5*VqLq*(nxJ.*FxSf1 + nyJ.*FySf1 + nzJ.*FzSf1));
rhs2 = 2*([VqPq VqLq]*r2 + .5*VqLq*(nxJ.*FxSf2 + nyJ.*FySf2 + nzJ.*FzSf2));
rhs3 = 2*([VqPq VqLq]*r3 + .5*VqLq*(nxJ.*FxSf3 + nyJ.*FySf3 + nzJ.*FzSf3));
rhs4 = 2*([VqPq VqLq]*r4 + .5*VqLq*(nxJ.*FxSf4 + nyJ.*FySf4 + nzJ.*FzSf4));
rhs5 = 2*([VqPq VqLq]*r5 + .5*VqLq*(nxJ.*FxSf5 + nyJ.*FySf5 + nzJ.*FzSf5));

disp('vol/surf averages')
[sum(sum(r1(1:Nq,:)))
sum(sum(r1(Nq+1:end,:)))]

% tmp = r1(1:Nq,:) + VqLq*(r1(Nq+1:end,:) + .5*(nxJ.*FxSf1 + nyJ.*FySf1 + nzJ.*FzSf1));
% sum(sum(tmp))
% sum(sum(VqLq*(.5*(nxJ.*FxSf1 + nyJ.*FySf1 + nzJ.*FzSf1))))

% wJq = diag(wq)*(Vq*J);
% 
% sum(sum(diag(wq)*(v1.*rhs1)))
% % compute entropy RHS - J already built in!
sum(sum(diag(wq)*(v1.*rhs1 + v2.*rhs2 + v3.*rhs3 + v4.*rhs4 + v5.*rhs5)))

%%

function [FxS1 FxS2 FxS3 FxS4 FxS5 FyS1 FyS2 FyS3 FyS4 FyS5 FzS1 FzS2 FzS3 FzS4 FzS5] = fluxes(rhoL,rhoR,uL,uR,vL,vR,wL,wR,EL,ER)

global gamma
global avg pfun betafun
betaL = betafun(rhoL,uL,vL,wL,EL);
betaR = betafun(rhoR,uR,vR,wR,ER);
rholog = logmean(rhoL, rhoR);
rhoavg = avg(rhoL, rhoR);
uavg = avg(uL, uR);
vavg = avg(vL, vR);
wavg = avg(wL, wR);
vnavg = 2.0 * (uavg .* uavg + vavg .* vavg + wavg .* wavg) - (avg(uL .* uL, uR .* uR) + avg(vL .* vL, vR .* vR) + avg(wL .* wL, wR .* wR));
beta_avg = avg(betaL, betaR);
pa = rhoavg ./ (2.0 * beta_avg);
f4aux = rholog ./ (2.0 * (gamma - 1.0) * logmean(betaL, betaR)) + pa + .5 * rholog .* vnavg;

FxS1 = rholog.*uavg;
FyS1 = rholog.*vavg;
FzS1 = rholog.*wavg;

FxS2 = FxS1.*uavg + pa;
FyS2 = FyS1.*uavg;
FzS2 = FzS1.*uavg;

FxS3 = FxS1.*vavg;
FyS3 = FyS1.*vavg + pa;
FzS3 = FzS1.*vavg;

FxS4 = FxS1.*wavg;
FyS4 = FyS1.*wavg;
FzS4 = FzS1.*wavg + pa;

FxS5 = f4aux.*uavg;
FyS5 = f4aux.*vavg;
FzS5 = f4aux.*wavg;

end
