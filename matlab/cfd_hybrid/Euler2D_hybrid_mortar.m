
% L2err_GLL_GLL = {};
% L2err_GQ_GLL = {};
% L2err_GQ_GQ = {};
% 
% for K1D = [6 12]
%     for N = 1:4
%         close all
%         
%         fprintf('N = %d, K1D = %d, GLL-GLL ------------\n',N,K1D)
%         
%         [rq1D_vol wq1D_vol] = JacobiGL(0,0,N);
%         [rq1D_face wq1D_face] = JacobiGL(0,0,N);
%         L2err = Euler2D(N, K1D, rq1D_vol, wq1D_vol, rq1D_face, wq1D_face)
%         L2err_GLL_GLL{N,K1D} = L2err;
%         
%         close all
%         
%         fprintf('N = %d, K1D = %d, GQ-GLL ------------\n',N,K1D)
%         
%         [rq1D_vol wq1D_vol] = JacobiGL(0,0,N);
%         [rq1D_face wq1D_face] = JacobiGQ(0,0,N);
%         L2err = Euler2D(N, K1D, rq1D_vol, wq1D_vol, rq1D_face, wq1D_face)
%         L2err_GLL_GQ{N,K1D} = L2err;
%         close all
%         
%         fprintf('N = %d, K1D = %d, GQ-GQ ------------\n',N,K1D)
%         
%         [rq1D_vol wq1D_vol] = JacobiGQ(0,0,N);
%         [rq1D_face wq1D_face] = JacobiGQ(0,0,N);
%         L2err = Euler2D(N, K1D, rq1D_vol, wq1D_vol, rq1D_face, wq1D_face)
%         L2err_GQ_GQ{N,K1D} = L2err;
%         
%     end
% end



function L2err = Euler2D(Nin, K1D, rq1D_quad_vol, wq1D_quad_vol, rq1D_face, wq1D_face)

useQuads = 0; mypath;
clear -globals
Globals2D

if nargin==0
    N = 3;
    K1D = 12; % K1D = multiples of 3
    
    [rq1D_quad_vol wq1D_quad_vol] = JacobiGQ(0,0,N);
    [rq1D_face wq1D_face] = JacobiGQ(0,0,N);
else
    N = Nin;
end

a = .125; %/K1D; % warping
CFL = .75;

plotMesh = 0;

FinalTime = 1.0;
global tau
tau = 1;

%% init hybrid mesh stuff

wadgProjEntropyVars = abs(a)>1e-8;

global QUAD TRI
QUAD = 1;
TRI = 2;

global x y xf yf xq yq
x = {}; y = {};
xf = {}; yf = {};
xq = {}; yq = {};

global DNr_type DNs_type VqPN_type VqLq_type
DNr_type = {}; DNs_type = {}; VqPN_type = {}; VqLq_type = {};
BNr = {}; BNs = {};

global mapP
mapP = {};

global rxJN_type sxJN_type ryJN_type syJN_type Jq_type nxJ_type nyJ_type sJ_type
rxJN_type = {}; sxJN_type = {};
ryJN_type = {}; syJN_type = {};
Jq_type = {};
sJ_type = {};

rxJ = {}; sxJ = {};
ryJ = {}; syJ = {};
nxJ_type = {}; nyJ_type = {};
J = {};

% for computing geofacs
V = {}; Vq = {}; Vf = {};
Dr = {}; Ds = {};
Vp = {}; xp = {}; yp = {};
nrJ = {}; nsJ = {};

%% setup quad mesh and tri mesh

[Nv, VX, VY, K, EToV] = QuadMesh2D(K1D,round(4/3*K1D));
iids = find(abs(abs(VX)-1)>1e-8 & abs(abs(VY)-1)>1e-8);
% VX(iids) = VX(iids) + a*randn(size(VX(iids)));
% VY(iids) = VY(iids) + a*randn(size(VX(iids)));
VX = VX + a*cos(pi/2*VX).*sin(pi*VY);
VY = VY + a*sin(pi*VX).*cos(pi/2*VY);

L = 7.5;
VX = (VX+1)/2 * L;
VY = VY*5;

% map nodes
r1 = [-1 1 1 -1]';
s1 = [-1 -1 1 1]';
r1D = JacobiGL(0,0,N);
[r s] = meshgrid(r1D); r = r(:); s = s(:);
V1 = Vandermonde2DQuad(1,r,s)/Vandermonde2DQuad(1,r1,s1);
x{QUAD} = V1*VX(EToV)';
y{QUAD} = V1*VY(EToV)';

% plot(x{QUAD},y{QUAD},'o')

% volume nodes
rq1D = rq1D_quad_vol;
wq1D = wq1D_quad_vol;
[rq sq] = meshgrid(rq1D); rq = rq(:); sq = sq(:);
[wrq wsq] = meshgrid(wq1D);
wq{QUAD} = wrq(:).*wsq(:);

rp1D = linspace(-1,1,10);
[rp sp] = meshgrid(rp1D); rp = rp(:); sp = sp(:);
V1D = Vandermonde1D(N,r1D);
D1D = GradVandermonde1D(N,r1D)/V1D;
V{QUAD} = Vandermonde2DQuad(N,r,s);

Vq{QUAD} = Vandermonde2DQuad(N,rq,sq)/V{QUAD};
xq{QUAD} = Vq{QUAD}*x{QUAD};
yq{QUAD} = Vq{QUAD}*y{QUAD};

Vp1D = Vandermonde1D(N,rp1D)/V1D;

Dr{QUAD} = kron(D1D,eye(N+1));
Ds{QUAD} = kron(eye(N+1),D1D);
Vp{QUAD} = kron(Vp1D,Vp1D);
xp{QUAD} = Vp{QUAD}*x{QUAD};
yp{QUAD} = Vp{QUAD}*y{QUAD};

% face nodes
rq1D = rq1D_face;
wq1D = wq1D_face;

e = ones(size(rq1D));
rfq = [rq1D; e; rq1D; -e];
sfq = [-e; rq1D; e; rq1D];
wfq = [wq1D; wq1D; wq1D; wq1D];
V1D = Vandermonde1D(N,JacobiGL(0,0,N));
Vfq = Vandermonde2DQuad(N,rfq,sfq)/V{QUAD};
Vf{QUAD} = Vfq;

nrJ{QUAD} = [0*e; e; 0*e; -e]; % sJ_type = 2 for all faces,
nsJ{QUAD} = [-e; 0*e; e; 0*e];

% quadrature operators
M = Vq{QUAD}'*diag(wq{QUAD})*Vq{QUAD};
Pq = M\(Vq{QUAD}'*diag(wq{QUAD}));
Lq = M\(Vfq'*diag(wfq));
VfPq{QUAD} = (Vfq*Pq);
VqPq{QUAD} = Vq{QUAD}*Pq;
VqLq = Vq{QUAD}*Lq;
VqLq_type{QUAD} = VqLq;
VqPN_type{QUAD} = [VqPq{QUAD} VqLq];

WN = diag([wq{QUAD};wfq]);

Qr = diag(wq{QUAD})*Vq{QUAD}*Dr{QUAD}*Pq;
Qs = diag(wq{QUAD})*Vq{QUAD}*Ds{QUAD}*Pq;
QNr = .5*[Qr-Qr' VfPq{QUAD}'*diag(wfq.*nrJ{QUAD});
    -diag(wfq.*nrJ{QUAD})*VfPq{QUAD} diag(0*wfq)];
QNs = .5*[Qs-Qs' VfPq{QUAD}'*diag(wfq.*nsJ{QUAD});
    -diag(wfq.*nsJ{QUAD})*VfPq{QUAD} diag(0*wfq)];

QNr(abs(QNr)<1e-8) = 0; QNr = sparse(QNr);
QNs(abs(QNs)<1e-8) = 0; QNs = sparse(QNs);

% make skew symmetric diff matrices
DNr_type{QUAD} = diag(1./[wq{QUAD};wfq])*QNr;
DNs_type{QUAD} = diag(1./[wq{QUAD};wfq])*QNs;

BNr{QUAD} = diag([zeros(size(Qr,1),1); nrJ{QUAD}]);
BNs{QUAD} = diag([zeros(size(Qr,1),1); nsJ{QUAD}]);

% make maps
xf{QUAD} = Vfq*x{QUAD};
yf{QUAD} = Vfq*y{QUAD};
[EToE,EToF]= tiConnectQuad2D(EToV);

[mapMq mapPq] = buildMaps(EToE,EToF,xf{QUAD},yf{QUAD});
mapP{QUAD} = mapPq;

if plotMesh
    for e = 1:K
        v = EToV(e,:);
        %         text(mean(VX(v)),mean(VY(v)),num2str(e))
        hold on
        v = [v v(1)];
        plot(VX(v),VY(v),'k-','linewidth',2)
    end
    axis equal
end
% return

%% build tri stuff

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D,round(4/3*K1D));
iids = find(abs(abs(VX)-1)>1e-8 & abs(abs(VY)-1)>1e-8);
% VX(iids) = VX(iids) + a*randn(size(VX(iids)));
% VY(iids) = VY(iids) + a*randn(size(VX(iids)));
VX = VX + a*cos(pi/2*VX).*sin(pi*VY);
VY = VY + a*sin(pi*VX).*cos(pi/2*VY);

VX = (VX+1)/2*L + L;
VY = VY*5;

[r s] = Nodes2D(N); [r s] = xytors(r,s);
r1 = [-1 1 -1]';
s1 = [-1 -1 1]';
V1 = Vandermonde2D(1,r,s)/Vandermonde2D(1,r1,s1);
x{TRI} = V1*VX(EToV)';
y{TRI} = V1*VY(EToV)';

V{TRI} = Vandermonde2D(N,r,s);
[Vr Vs] = GradVandermonde2D(N,r,s);
Dr{TRI} = Vr/V{TRI};
Ds{TRI} = Vs/V{TRI};

[rq sq wqtri] = Cubature2D(2*N);
wq{TRI} = wqtri;
Vq{TRI} = Vandermonde2D(N,rq,sq)/V{TRI};
xq{TRI} = Vq{TRI}*x{TRI};
yq{TRI} = Vq{TRI}*y{TRI};

[rp sp] = EquiNodes2D(5); [rp sp] = xytors(rp,sp);
Vp{TRI} = Vandermonde2D(N,rp,sp)/V{TRI};
xp{TRI} = Vp{TRI}*x{TRI};
yp{TRI} = Vp{TRI}*y{TRI};

% face nodes
rq1D = rq1D_face;
wq1D = wq1D_face;

rfq = [rq1D; -rq1D; -ones(size(rq1D))];
sfq = [-ones(size(rq1D)); rq1D; -rq1D];
wfq = [wq1D;wq1D;wq1D];
Vfq = Vandermonde2D(N,rfq,sfq)/V{TRI};
Vf{TRI} = Vfq;
nrJ{TRI} = [-zeros(size(rq1D)); ones(size(rq1D)); -ones(size(rq1D)); ];
nsJ{TRI} = [-ones(size(rq1D)); ones(size(rq1D)); -zeros(size(rq1D)); ];

% operators
M = Vq{TRI}'*diag(wq{TRI})*Vq{TRI};
Pq = M\(Vq{TRI}'*diag(wq{TRI}));
Lq = M\(Vfq'*diag(wfq));
VfPq{TRI} = (Vfq*Pq);
VqPq{TRI} = Vq{TRI}*Pq;
VqLq = Vq{TRI}*Lq;
VqLq_type{TRI} = VqLq;
VqPN_type{TRI} = [VqPq{TRI} VqLq];

WN = diag([wq{TRI};wfq]);
Qr = diag(wq{TRI})*Vq{TRI}*Dr{TRI}*Pq;
Qs = diag(wq{TRI})*Vq{TRI}*Ds{TRI}*Pq;
QNr = .5*[Qr-Qr' VfPq{TRI}'*diag(wfq.*nrJ{TRI});
    -diag(wfq.*nrJ{TRI})*VfPq{TRI} diag(0*wfq)];
QNs = .5*[Qs-Qs' VfPq{TRI}'*diag(wfq.*nsJ{TRI});
    -diag(wfq.*nsJ{TRI})*VfPq{TRI} diag(0*wfq)];
QNr(abs(QNr)<1e-8) = 0; QNr = sparse(QNr);
QNs(abs(QNs)<1e-8) = 0; QNs = sparse(QNs);

% make skew symmetric diff matrices
DNr_type{TRI} = diag(1./[wq{TRI};wfq])*QNr;
DNs_type{TRI} = diag(1./[wq{TRI};wfq])*QNs;

BNr{TRI} = diag([zeros(size(Qr,1),1); nrJ{TRI}]);
BNs{TRI} = diag([zeros(size(Qr,1),1); nsJ{TRI}]);

% make maps
xf{TRI} = Vfq*x{TRI};
yf{TRI} = Vfq*y{TRI};
[EToE,EToF]= tiConnect2D(EToV);
[mapMq mapPq] = buildMaps(EToE,EToF,xf{TRI},yf{TRI});
mapP{TRI} = mapPq;

if plotMesh
    hold on
    for e = 1:K
        v = EToV(e,:);
        v = [v v(1)];
        plot(VX(v),VY(v),'k-','linewidth',2)
        hold on
    end
    axis off
    return
end


%% build geometric factors

wJq_type = {};
for TYPE = 1:2
    
    [rx,sx,ry,sy,Jt] = GeometricFactors2D(x{TYPE},y{TYPE},Dr{TYPE},Ds{TYPE});
    rxJ{TYPE} = rx.*Jt;
    sxJ{TYPE} = sx.*Jt;
    ryJ{TYPE} = ry.*Jt;
    syJ{TYPE} = sy.*Jt;
    J{TYPE} = Jt;
    
    % used in computation
    VN = [Vq{TYPE}; Vf{TYPE}];
    rxJN_type{TYPE} = VN*rxJ{TYPE};    sxJN_type{TYPE} = VN*sxJ{TYPE};
    ryJN_type{TYPE} = VN*ryJ{TYPE};    syJN_type{TYPE} = VN*syJ{TYPE};
    Jq_type{TYPE} = Vq{TYPE}*J{TYPE};
    
    dnrJ = spdiag(nrJ{TYPE});
    dnsJ_type = spdiag(nsJ{TYPE});
    nxJ_type{TYPE} = dnrJ*(Vf{TYPE}*rxJ{TYPE}) + dnsJ_type*(Vf{TYPE}*sxJ{TYPE});
    nyJ_type{TYPE} = dnrJ*(Vf{TYPE}*ryJ{TYPE}) + dnsJ_type*(Vf{TYPE}*syJ{TYPE});
    sJ_type{TYPE} = sqrt(nxJ_type{TYPE}.^2 + nyJ_type{TYPE}.^2);
    
    wJq_type{TYPE} = diag(wq{TYPE})*(Vq{TYPE}*Jt);
    
    %     quiver(xf{TYPE},yf{TYPE},nxJ_type{TYPE},nyJ_type{TYPE})
    %     hold on
    %     axis equal
    
end

%% test differentiation
%
% if 0
%     f = @(x,y) exp(x+y);
%     df = @(x,y) exp(x+y);
%
%     for TYPE = 1:2
%         u = f(x{TYPE},y{TYPE});
%         dudxJ = rxJ{TYPE}.*(Dr{TYPE}*u) + sxJ{TYPE}.*(Ds{TYPE}*u);
%         norm(dudxJ - df(x{TYPE},y{TYPE}).*J{TYPE},'fro')
%
%         uN = [Vq{TYPE};Vf{TYPE}]*u;
%         dudxJ = VqPN_type{TYPE}*(rxJN_type{TYPE}.*((.5*BNr{TYPE} + DNr_type{TYPE})*uN) + sxJN_type{TYPE}.*((.5*BNs{TYPE} + DNs_type{TYPE})*uN));
%         norm(dudxJ - df(xq{TYPE},yq{TYPE}).*(Vq{TYPE}*J{TYPE}),'fro')
%     end
% end

%% make maps for quad/tri interface

xQmin = min(x{QUAD}(:)); xQmax = max(x{QUAD}(:));
yQmin = min(y{QUAD}(:)); yQmax = max(y{QUAD}(:));
xTmin = min(x{TRI}(:)); xTmax = max(x{TRI}(:));
yTmin = min(y{TRI}(:)); yTmax = max(y{TRI}(:));

global mapMT mapPT mapMQ mapPQ
Nfaces = 4;
K = size(xf{QUAD},2);
Nfp = size(xf{QUAD},1)/Nfaces;
xfQ = reshape(xf{QUAD},Nfp,Nfaces*K);
yfQ = reshape(yf{QUAD},Nfp,Nfaces*K);
mapMQ = reshape(1:Nfp*Nfaces*K,Nfp,Nfaces*K);

Nfaces = 3;
K = size(xf{TRI},2);
xfT = reshape(xf{TRI},Nfp,Nfaces*K);
yfT = reshape(yf{TRI},Nfp,Nfaces*K);
mapMT = reshape(1:Nfp*Nfaces*K,Nfp,Nfaces*K);

% make quad/tri maps periodic in y-direction
fQ1 = find(sum(abs(yfQ-yQmax)) < 1e-9); % top
fQ2 = find(sum(abs(yfQ-yQmin)) < 1e-9); % bottom
mapP{QUAD}(mapMQ(:,fQ1)) = mapMQ(:,fQ2); % ASSUMES SAME ORDERING
mapP{QUAD}(mapMQ(:,fQ2)) = mapMQ(:,fQ1);

fT1 = find(sum(abs(yfT-yTmax)) < 1e-9); % top
fT2 = find(sum(abs(yfT-yTmin)) < 1e-9); % bottom
mapP{TRI}(mapMT(:,fT1)) = flipud(mapMT(:,fT2)); % ASSUMES LOCALLY REVERSED ORDERING
mapP{TRI}(mapMT(:,fT2)) = flipud(mapMT(:,fT1));

% make quad/tri maps periodic in x-direction (TESTING ONLY)
fQ1 = find(sum(abs(xfQ-xQmin)) < 1e-9); % left
fQ2 = find(sum(abs(xfQ-xQmax)) < 1e-9); % right
mapP{QUAD}(mapMQ(:,fQ1)) = (mapMQ(:,fQ2)); % ASSUMES SAME ORDERING
mapP{QUAD}(mapMQ(:,fQ2)) = (mapMQ(:,fQ1));

fT1 = find(sum(abs(xfT-xTmin)) < 1e-9); % left
fT2 = find(sum(abs(xfT-xTmax)) < 1e-9); % right
mapP{TRI}(mapMT(:,fT1)) = flipud(mapMT(:,fT2)); % ASSUMES LOCALLY REVERSED ORDERING
mapP{TRI}(mapMT(:,fT2)) = flipud(mapMT(:,fT1));

% clf
% plot(xfT,yfT,'o')
% hold on
% text(xfT(:)-.025,yfT(:)-.025,num2str((1:length(xfT(:)))'))
% plot(xfQ,yfQ,'x')
% text(xfQ(:)+.025,yfQ(:)+.025,num2str((1:length(xfQ(:)))'))
% return

if 1
    % make hybrid tri-quad maps
    fQ1 = find(sum(abs(xfQ - xQmax)) < 1e-9); % dividing line
    fQ2 = find(sum(abs(xfQ - xQmin)) < 1e-9);
    fQ = [fQ1(:); fQ2(:)];
    mapMQ = mapMQ(:,fQ);
    
    fT1 = find(sum(abs(xfT - xTmin)) < 1e-9);
    fT2 = find(sum(abs(xfT - xTmax)) < 1e-9); % dividing line
    fT = [fT1(:); fT2(:)];
    mapMT = mapMT(:,fT);
    
    xfQ_hybrid = xfQ(mapMQ);
    yfQ_hybrid = yfQ(mapMQ);
    xfT_hybrid = xfT(mapMT);
    yfT_hybrid = yfT(mapMT);
    
    if norm(mean(yfQ_hybrid,1)-mean(yfT_hybrid,1))>1e-8
        % if gets here, need to check ordering!
        keyboard
    end
    
    mapPQ = zeros(Nfp,length(fQ));
    mapPT = zeros(Nfp,length(fT));
    for f = 1:length(fQ)
        [x1 y1] = meshgrid(xfQ_hybrid(:,f),yfQ_hybrid(:,f));
        [x2 y2] = meshgrid(xfT_hybrid(:,f),yfT_hybrid(:,f));
        %D = abs(x1-x2') + abs(y1-y2');
        D = abs(y1-y2'); % WARNING: ASSUMES ORDERED + CARTESIAN INTERFACES INCLUDING PERIODIC BOUNDARIES
        
        [p,~] = find(D < 1e-9);
        mapPQ(:,f) = p + (f-1)*Nfp;
        
        [p,~] = find(D' < 1e-9);
        mapPT(:,f) = p + (f-1)*Nfp;
    end
    
    % check matches
    u = {};
    for TYPE = 1:2
        u{TYPE} = 1+sin(pi/L*(x{TYPE} + y{TYPE}));
        uf{TYPE} = Vf{TYPE}*u{TYPE};
    end
    
    mapPT = mapMQ(mapPT);
    mapPQ = mapMT(mapPQ);
    err = norm(uf{TRI}(mapMT) - uf{QUAD}(mapPT),'fro') + norm(uf{QUAD}(mapMQ) - uf{TRI}(mapPQ),'fro');
    if err > 1e-12
        keyboard
    end
end

keyboard


%% visualize connectivity

% for i = 1:length(mapP{TRI}(:))
%     clf
%     hold on
%     plot(xf{TRI},yf{TRI},'o')
%     plot(xf{TRI}(i),yf{TRI}(i),'+','markersize',20)
%     plot(xf{TRI}(mapP{TRI}(i)),yf{TRI}(mapP{TRI}(i)),'x','markersize',20)
%     pause%(.1)
% end
%
% for i = 1:length(mapP{QUAD}(:))
%     clf
%     hold on
%     plot(xf{QUAD},yf{QUAD},'o')
%     plot(xf{QUAD}(i),yf{QUAD}(i),'+','linewidth',2,'markersize',20)
%     plot(xf{QUAD}(mapP{QUAD}(i)),yf{QUAD}(mapP{QUAD}(i)),'x','linewidth',2,'markersize',20)
%     pause%(.1)
% end
% return

% for i = 1:length(mapPQ(:))
%     clf
%     plot(xfT,yfT,'o')
%     hold on
%     plot(xfQ,yfQ,'x')
%     plot(xfT(mapMT(i)),yfT(mapMT(i)),'s','markersize',20)
%     plot(xfQ(mapPT(i)),yfQ(mapPT(i)),'*','markersize',20)
%     pause(.5)
% end
% return


%% define fluxes

global gamma
gamma = 1.4;

global U1 U2 U3 U4 V1 V2 V3 V4
rhoe = @(rho,rhou,rhov,E) E - .5*(rhou.^2+rhov.^2)./rho;
pcons = @(rho,rhou,rhov,E) (gamma-1)*rhoe(rho,rhou,rhov,E);
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


global pfun beta pavg plogmean vnormavg avg

avg = @(x,y) .5*(x+y);
pfun = @(rho,u,v,E) (gamma-1)*(E-.5*rho.*(u.^2+v.^2));
beta = @(rho,u,v,E) rho./(2*pfun(rho,u,v,E));
pavg     = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER)     avg(rhoL,rhoR)./(2*avg(beta(rhoL,uL,vL,EL),beta(rhoR,uR,vR,ER)));
plogmean = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR)./(2*logmean(beta(rhoL,uL,vL,EL),beta(rhoR,uR,vR,ER)));
vnormavg = @(uL,vL,uR,vR) 2*(avg(uL,uR).^2 + avg(vL,vR).^2) - (avg(uL.^2,uR.^2) + avg(vL.^2,vR.^2));


%% run solver

rk4a = [            0.0 ...
    -567301805773.0/1357537059087.0 ...
    -2404267990393.0/2016746695238.0 ...
    -3550918686646.0/2091501179385.0  ...
    -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
    5161836677717.0/13612068292357.0 ...
    1720146321549.0/2090206949498.0  ...
    3134564353537.0/4481467310338.0  ...
    2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
    1432997174477.0/9575080441755.0 ...
    2526269341429.0/6820363962896.0 ...
    2006345519317.0/3224310063776.0 ...
    2802321613138.0/2924317926251.0];


for TYPE = 1:2
    
    [rhoex uex vex pex] = vortexSolution(xq{TYPE},yq{TYPE},0);
    rho{TYPE} = VqPq{TYPE}*rhoex;
    rhou{TYPE} = VqPq{TYPE}*(rhoex.*uex);
    rhov{TYPE} = VqPq{TYPE}*(rhoex.*vex);
    E{TYPE} = VqPq{TYPE}*(pex/(gamma-1) + .5*rhoex.*(uex.^2+vex.^2));
    
    rhoq{TYPE} = VqPq{TYPE}*rhoex;
    uq{TYPE} = VqPq{TYPE}*(uex);
    vq{TYPE} = VqPq{TYPE}*(vex);
    Eq{TYPE} = VqPq{TYPE}*(pex/(gamma-1) + .5*rhoex.*(uex.^2+vex.^2));
    
    rhoM{TYPE} = [];
    uM{TYPE} = [];
    vM{TYPE} = [];
    EM{TYPE} = [];
    
    % Runge-Kutta residual storage
    K = size(xq{TYPE},2);
    Nq = size(xq{TYPE},1);
    res1{TYPE} = zeros(Nq,K);
    res2{TYPE} = zeros(Nq,K);
    res3{TYPE} = zeros(Nq,K);
    res4{TYPE} = zeros(Nq,K);
end

% [a b c d] = vortexSolution(0,0,0)
% clf
% for TYPE=1:2
%     vv = rho{TYPE};
%     color_line3(xq{TYPE},yq{TYPE},vv,vv,'.')
%     hold on
% end
% axis equal
% return

% compute time step size
h_estimate = min(min(J{TYPE})./max(sJ_type{TYPE})); % estimate h using ||J||/||sJ||
CN = (N+1)*(N+2)/2; % trace const for quads (larger than tris)
dt = CFL * h_estimate/CN; % h = J/Jf
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

entropy = zeros(Nsteps,1);
rhstest = zeros(Nsteps,1);
v1 = {};
v2 = {};
v3 = {};
v4 = {};
v1M = {};
v2M = {};
v3M = {};
v4M = {};
            
figure(1)
for i = 1:Nsteps
    for INTRK = 1:5
        
        for TYPE = 1:2
            q1q = VqPq{TYPE}*V1(rho{TYPE},rhou{TYPE},rhov{TYPE},E{TYPE});
            q2q = VqPq{TYPE}*V2(rho{TYPE},rhou{TYPE},rhov{TYPE},E{TYPE});
            q3q = VqPq{TYPE}*V3(rho{TYPE},rhou{TYPE},rhov{TYPE},E{TYPE});
            q4q = VqPq{TYPE}*V4(rho{TYPE},rhou{TYPE},rhov{TYPE},E{TYPE});
            
            v1{TYPE} = q1q;
            v2{TYPE} = q2q;
            v3{TYPE} = q3q;
            v4{TYPE} = q4q;
            
            v1M{TYPE} = VfPq{TYPE}*q1q;
            v2M{TYPE} = VfPq{TYPE}*q2q;
            v3M{TYPE} = VfPq{TYPE}*q3q;
            v4M{TYPE} = VfPq{TYPE}*q4q;
            
            rhoq{TYPE}  = U1(q1q,q2q,q3q,q4q);
            rhouq       = U2(q1q,q2q,q3q,q4q);
            rhovq       = U3(q1q,q2q,q3q,q4q);
            Eq{TYPE}    = U4(q1q,q2q,q3q,q4q);
            uq{TYPE} = rhouq./rhoq{TYPE};
            vq{TYPE} = rhovq./rhoq{TYPE};                        
        end
        
        for TYPE = 1:2
            
            [rhs1, rhs2, rhs3, rhs4]  = RHS2Dsimple(rhoq{TYPE},uq{TYPE},vq{TYPE},Eq{TYPE},...
                v1M,v2M,v3M,v4M,TYPE);
            
            if (INTRK==5)
                rhstest(i) = rhstest(i) + sum(sum(wJq_type{TYPE}.*(rhs1.*v1{TYPE} + rhs2.*v2{TYPE} + rhs3.*v3{TYPE} + rhs4.*v4{TYPE})));
            end
            
            res1{TYPE} = rk4a(INTRK)*res1{TYPE} + dt*rhs1;
            res2{TYPE} = rk4a(INTRK)*res2{TYPE} + dt*rhs2;
            res3{TYPE} = rk4a(INTRK)*res3{TYPE} + dt*rhs3;
            res4{TYPE} = rk4a(INTRK)*res4{TYPE} + dt*rhs4;
            
            rho{TYPE}  = rho{TYPE}  + rk4b(INTRK)*res1{TYPE};
            rhou{TYPE} = rhou{TYPE} + rk4b(INTRK)*res2{TYPE};
            rhov{TYPE} = rhov{TYPE} + rk4b(INTRK)*res3{TYPE};
            E{TYPE}    = E{TYPE}    + rk4b(INTRK)*res4{TYPE};
        end
    end;
    
    for TYPE = 1:2
        Sq = -rho{TYPE}.*s(rho{TYPE},rhou{TYPE},rhov{TYPE},E{TYPE});
        entropy(i) = entropy(i) + sum(sum(wJq_type{TYPE}.*Sq));
    end
    
    if  (mod(i,10)==0 || i==Nsteps)
        clf
        hold on
        for TYPE=1:2
            vv = rho{TYPE};
            color_line3(xq{TYPE},yq{TYPE},vv,vv,'.');
        end
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f, step %d out of %d, N = %d, K1D = %d',dt*i,i,Nsteps,N,K1D))
        drawnow
    end
end

figure
entropy = entropy - min(entropy) + 1;
semilogy(dt*(1:Nsteps),abs(rhstest),'o--')
hold on
semilogy(dt*(1:Nsteps),entropy,'x--')
legend('Entropy rhs','entropy')

L2err = 0;
for TYPE=1:2
    if TYPE==1 %quad
        [rq1D2 wq1D2] = JacobiGQ(0,0,N+2);
        [rq2 sq2] = meshgrid(rq1D2); rq2 = rq2(:); sq2 = sq2(:);
        [wrq wsq] = meshgrid(wq1D2); wq2 = wrq(:).*wsq(:);
        Vq2 = Vandermonde2DQuad(N,rq2,sq2)/V{QUAD};
    else
        [rq2 sq2 wq2] = Cubature2D(2*N+2);
        Vq2 = Vandermonde2D(N,rq2,sq2)/V{TRI};
    end
    VqPq = Vq2*((Vq{TYPE}'*diag(wq{TYPE})*Vq{TYPE})\(Vq{TYPE}'*diag(wq{TYPE})));
    xq2 = Vq2*x{TYPE};
    yq2 = Vq2*y{TYPE};
    wJq2 = diag(wq2)*(Vq2*J{TYPE});
    [rhoex uex vex pex] = vortexSolution(xq2,yq2,FinalTime);
    rhouex = rhoex.*uex;
    rhovex = rhoex.*vex;
    Eex = pex/(gamma-1) + .5*rhoex.*(uex.^2+vex.^2);
    
    err = wJq2.*((VqPq*rho{TYPE}-rhoex).^2 + (VqPq*rhou{TYPE}-rhouex).^2 + (VqPq*rhov{TYPE}-rhovex).^2 + (VqPq*E{TYPE}-Eex).^2);
    L2err = L2err + sum(err(:));
end
L2err = sqrt(L2err);

end
%%

function [rhs1 rhs2 rhs3 rhs4] = RHS2Dsimple(rhoq,uq,vq,Eq,v1_type,v2_type,v3_type,v4_type,TYPE) %rhoM_type,uM_type,vM_type,EM_type,TYPE)

% Globals2D;

global DNr_type DNs_type VqPN_type VqLq_type
global mapP mapMT mapPT mapMQ mapPQ
global rxJN_type sxJN_type ryJN_type syJN_type Jq_type nxJ_type nyJ_type sJ_type
global beta avg

DNr = DNr_type{TYPE};
DNs = DNs_type{TYPE};
VqPN = VqPN_type{TYPE};
VqLq = VqLq_type{TYPE};

rxJN = rxJN_type{TYPE};
sxJN = sxJN_type{TYPE};
ryJN = ryJN_type{TYPE};
syJN = syJN_type{TYPE};
Jq = Jq_type{TYPE};
nxJ = nxJ_type{TYPE};
nyJ = nyJ_type{TYPE};
sJ = sJ_type{TYPE};

K = size(Jq,2); % cols

v1M = v1_type{TYPE};
v2M = v2_type{TYPE};
v3M = v3_type{TYPE};
v4M = v4_type{TYPE};

v1P = v1M(mapP{TYPE});
v2P = v2M(mapP{TYPE});
v3P = v3M(mapP{TYPE});
v4P = v4M(mapP{TYPE});

% hybrid coupling
global QUAD TRI
if 1
    % get nbr info
    if TYPE==QUAD
        v1P(mapMQ) = v1_type{TRI}(mapPQ);
        v2P(mapMQ) = v2_type{TRI}(mapPQ);
        v3P(mapMQ) = v3_type{TRI}(mapPQ);
        v4P(mapMQ) = v4_type{TRI}(mapPQ);
    elseif TYPE==TRI
        v1P(mapMT) = v1_type{QUAD}(mapPT);
        v2P(mapMT) = v2_type{QUAD}(mapPT);
        v3P(mapMT) = v3_type{QUAD}(mapPT);
        v4P(mapMT) = v4_type{QUAD}(mapPT);
    end
end

global U1 U2 U3 U4 V1 V2 V3 V4
rhoM  = U1(v1M,v2M,v3M,v4M);
rhouM = U2(v1M,v2M,v3M,v4M);
rhovM = U3(v1M,v2M,v3M,v4M);
EM    = U4(v1M,v2M,v3M,v4M);
uM = rhouM./rhoM;
vM = rhovM./rhoM;

rhoP  = U1(v1P,v2P,v3P,v4P);
rhouP = U2(v1P,v2P,v3P,v4P);
rhovP = U3(v1P,v2P,v3P,v4P);
EP    = U4(v1P,v2P,v3P,v4P);
uP = rhouP./rhoP;
vP = rhovP./rhoP;

betaq = beta(rhoq,uq,vq,Eq);
betaM = beta(rhoM,uM,vM,EM);
betaP = beta(rhoP,uP,vP,EP);

% Lax-Friedrichs flux
global gamma tau
unorm2 = (uM.^2+vM.^2);
pM = (gamma-1)*(EM - .5*rhoM.*unorm2);
cvel = sqrt(gamma*pM./rhoM);
lam = sqrt(unorm2)+cvel;

unorm2 = (uP.^2+vP.^2);
pP = (gamma-1)*(EP - .5*rhoP.*unorm2);
cvel = sqrt(gamma*pP./rhoP);
lamP = sqrt(unorm2)+cvel;

LFc = max(lamP,lam);
dQ1 = rhoP-rhoM;
dQ2 = rhoP.*uP-rhoM.*uM;
dQ3 = rhoP.*vP-rhoM.*vM;
dQ4 = EP-EM;
Lf1 = tau*LFc.*dQ1.*sJ;
Lf2 = tau*LFc.*dQ2.*sJ;
Lf3 = tau*LFc.*dQ3.*sJ;
Lf4 = tau*LFc.*dQ4.*sJ;

[FxS1 FxS2 FxS3 FxS4 FyS1 FyS2 FyS3 FyS4] = euler_flux(rhoM,rhoP,uM,uP,vM,vP,betaM,betaP,gamma);

fSf1 = nxJ.*FxS1 + nyJ.*FyS1;
fSf2 = nxJ.*FxS2 + nyJ.*FyS2;
fSf3 = nxJ.*FxS3 + nyJ.*FyS3;
fSf4 = nxJ.*FxS4 + nyJ.*FyS4;

fSf1 = fSf1  - .25*Lf1;
fSf2 = fSf2  - .25*Lf2;
fSf3 = fSf3  - .25*Lf3;
fSf4 = fSf4  - .25*Lf4;

%% mortar 


% for ff = 1:length(split_faces)
%     f = split_faces(ff);
%     ncids = Nfaces*K + [2*ff-1 2*ff];
%     
%     % flux vs mortar points
%     rhom = reshape(rhoM(:,ncids),2*Nfp,1);
%     um = reshape(uM(:,ncids),2*Nfp,1);
%     vm = reshape(vM(:,ncids),2*Nfp,1);
%     betam = reshape(betaM(:,ncids),2*Nfp,1);
%     Em = reshape(EM(:,ncids),2*Nfp,1);
%     
%     [rhox, rhoy] = meshgrid(rhom,rhoM(:,f));
%     [ux, uy] = meshgrid(um,uM(:,f));
%     [vx, vy] = meshgrid(vm,vM(:,f));
%     [betax, betay] = meshgrid(betam,betaM(:,f));
%     [Ex, Ey] = meshgrid(Em,EM(:,f));
%     
%     [nxJ1, nxJ2] = meshgrid(reshape(nxJ(:,ncids),2*Nfp,1),nxJ(:,f));
%     [nyJ1, nyJ2] = meshgrid(reshape(nyJ(:,ncids),2*Nfp,1),nyJ(:,f));
%     nxJa = .5*(nxJ1+nxJ2);
%     nyJa = .5*(nyJ1+nyJ2);
%     
%     % flux correction terms: E
%     [FxS1 FxS2 FxS3 FxS4 FyS1 FyS2 FyS3 FyS4] = euler_flux(rhox,rhoy,ux,uy,vx,vy,Ex,Ey,betax,betay,gamma);
%     
%     fSc1 = nxJa.*FxS1 + nyJa.*FyS1;
%     fSc2 = nxJa.*FxS2 + nyJa.*FyS2;
%     fSc3 = nxJa.*FxS3 + nyJa.*FyS3;
%     fSc4 = nxJa.*FxS4 + nyJa.*FyS4;
%     
%     fc1 = sum(Emf.*fSc1,2) - Emf*sum(Efm.*fSc1',2);
%     fc2 = sum(Emf.*fSc2,2) - Emf*sum(Efm.*fSc2',2);
%     fc3 = sum(Emf.*fSc3,2) - Emf*sum(Efm.*fSc3',2);
%     fc4 = sum(Emf.*fSc4,2) - Emf*sum(Efm.*fSc4',2);
%     
%     fSf1(:,f) = Emf*reshape(fSf1(:,ncids),2*Nfp,1) + fc1;
%     fSf2(:,f) = Emf*reshape(fSf2(:,ncids),2*Nfp,1) + fc2;
%     fSf3(:,f) = Emf*reshape(fSf3(:,ncids),2*Nfp,1) + fc3;
%     fSf4(:,f) = Emf*reshape(fSf4(:,ncids),2*Nfp,1) + fc4;
% end

%%

rhoN = [rhoq; rhoM];
uN = [uq; uM];
vN = [vq; vM];
betaN = [betaq;betaM];

divF1 = zeros(size(DNr,1),K);
divF2 = zeros(size(DNr,1),K);
divF3 = zeros(size(DNr,1),K);
divF4 = zeros(size(DNr,1),K);
for e = 1:K
    [rhox, rhoy] = meshgrid(rhoN(:,e));
    [ux, uy] = meshgrid(uN(:,e));
    [vx, vy] = meshgrid(vN(:,e));
    [betax, betay] = meshgrid(betaN(:,e));
    
    % optimized evaluations
    [FxS1 FxS2 FxS3 FxS4 FyS1 FyS2 FyS3 FyS4] = euler_flux(rhox,rhoy,ux,uy,vx,vy,betax,betay,gamma);
    
    [rxJ1, rxJ2] = meshgrid(rxJN(:,e));
    [sxJ1, sxJ2] = meshgrid(sxJN(:,e));
    [ryJ1, ryJ2] = meshgrid(ryJN(:,e));
    [syJ1, syJ2] = meshgrid(syJN(:,e));
    rxJK = avg(rxJ1,rxJ2);  sxJK = avg(sxJ1,sxJ2);
    ryJK = avg(ryJ1,ryJ2);  syJK = avg(syJ1,syJ2);
    
    Dx = DNr.*rxJK + DNs.*sxJK;
    Dy = DNr.*ryJK + DNs.*syJK;
    
    divF1(:,e) = sum(Dx.*FxS1,2) + sum(Dy.*FyS1,2);
    divF2(:,e) = sum(Dx.*FxS2,2) + sum(Dy.*FyS2,2);
    divF3(:,e) = sum(Dx.*FxS3,2) + sum(Dy.*FyS3,2);
    divF4(:,e) = sum(Dx.*FxS4,2) + sum(Dy.*FyS4,2);
    
end

rhs1 = 2*VqPN*divF1 + VqLq*(fSf1);
rhs2 = 2*VqPN*divF2 + VqLq*(fSf2);
rhs3 = 2*VqPN*divF3 + VqLq*(fSf3);
rhs4 = 2*VqPN*divF4 + VqLq*(fSf4);

% assume mass lumping or Jq const (tris)
rhs1 = -(rhs1./Jq);
rhs2 = -(rhs2./Jq);
rhs3 = -(rhs3./Jq);
rhs4 = -(rhs4./Jq);

end


function [FxS1 FxS2 FxS3 FxS4 FyS1 FyS2 FyS3 FyS4] = euler_flux(rhox,rhoy,ux,uy,vx,vy,betax,betay,gamma)

rholog = logmean(rhox,rhoy);

% arithmetic avgs
rhoavg = .5*(rhox+rhoy);
uavg = .5*(ux+uy);
vavg = .5*(vx+vy);

vnavg = 2*(uavg.^2 + vavg.^2) - .5*((ux.^2+uy.^2) + (vx.^2+vy.^2));
pa = rhoavg./(betax+betay);

FxS1 = rholog.*uavg;      FyS1 = rholog.*vavg;
FxS2 = FxS1.*uavg + pa;   FyS2 = FyS1.*uavg;
FxS3 = FyS2;              FyS3 = FyS1.*vavg + pa;
f4aux = rholog./(2*(gamma-1)*logmean(betax,betay)) + pa + .5*rholog.*vnavg;
FxS4 = f4aux.*uavg;
FyS4 = f4aux.*vavg;

end

function [rho u v p] = vortexSolution(x,y,t)

% [0 20] x [-5 5]
% x = (1+x) * 10;
% y = y * 10;
% t = t * 10;

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

if 0
    % constant free stream test
    rho = 2 + 0*x;
    u = .2*rho;
    v = .1*rho;
    p = rho.^gamma;
end

if 0
    % pulse condition
    x0 = 2.5;
    rho = 2 + (abs(x-x0) < 5);
    u = 0*rho;
    v = 0*rho;
    p = rho.^gamma;
    
end

end

function [mapMq mapPq] = buildMaps(EToE,EToF,xf,yf)

NfqNfaces = size(xf,1);
Nfaces = size(EToF,2);
Nfq = NfqNfaces/Nfaces;
K = size(EToE,1);
mapMq = reshape(1:length(xf(:)),Nfq*Nfaces,K);
mapPq = mapMq;

tol = 1e-11;
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
            
            %             % NOTE - does not work if K1D is too small!!
            %             if length(p) == 0
            %                 % assume periodic boundary, find match in x,y
            %                 [px,~] = find(DX<tol);
            %                 [py,~] = find(DY<tol);
            %                 if length(px)==0
            %                     p = py;
            %                 elseif length(py)==0
            %                     p = px;
            %                 else
            %                     keyboard
            %                 end
            %             end
            fids = id2(p) + (enbr-1)*(Nfq*Nfaces);
            mapPq(id1,e) = fids(:);
        end
    end
end

end
