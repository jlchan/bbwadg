clear
Globals2D

%% estimate correct K

% % 6144 for lowest order
% 
% N = 1;
% K1D = 2.^(1:6)
% 
% K = 2*K1D.^2;
% Np = (N+1)*(N+2)/2;
% NpK = Np*K;
% K1D(find(NpK > 6144,1,'first'))


%%

a = .05; % warping factor
% a = 1/8;
% FinalTime = 10;

errAffine = [0.0040 7.1861e-04 1.7965e-04 3.8038e-05 6.2680e-06 3.4493e-06];
errCurvedRef = [0.0040  0.0011 7.9894e-04 5.5982e-04 3.5219e-04 3.5968e-04];
errCurvedPhys = [0.0040 7.4954e-04 2.2185e-04 4.6794e-05 1.0430e-05 6.1943e-06];
% errCurvedPhys2Np2 = [0.0040  7.4954e-04 2.4793e-04 ];
semilogy(errCurvedRef(1:6),'bs--','linewidth',2,'markersize',12,'MarkerFaceColor',[.49 1 .63],'DisplayName','DG, coarsened curved mesh')
hold on
semilogy(errAffine(1:6),'ko--','linewidth',2,'markersize',12,'MarkerFaceColor',[.49 1 .63],'DisplayName','DG, coarsened affine mesh')
semilogy(errCurvedPhys(1:6),'r^--','linewidth',2,'markersize',12,'MarkerFaceColor',[.49 1 .63],'DisplayName','Phys.\ frame DG, curved mesh')
grid on
fontsize = 20;
set(gca,'fontsize',fontsize,'TickLabelInterpreter','latex');
xlabel('Polynomial degree $N$','fontsize',fontsize,'interpreter','latex')
ylabel('$L^2$ error','fontsize',fontsize,'interpreter','latex')
h = legend('show','location','southwest');
set(h,'interpreter','latex')
print(gcf,'-dpng','~/Desktop/proposals/CAREER/figs/physframe_convergence.png')
% export_fig ~/Desktop/proposals/CAREER/figs/physframe_convergence.png -transparent -painters
return

% N =   1,   2,  3,  4,  5, 6
% K1D = 32, 22, 17, 14, 12, 10

N = 6;
K1D = 10;
FinalTime = 1;
plotMesh = 1;
usePhysOps = 1;
Ngeo = 1;

wadgProjEntropyVars = abs(a)>1e-8;
CFL = 1.75;
global tau
tau = 1;

Lx = 1; Ly = 1; ratiox = 1; ratioy = 1;
% Lx = 7.5; Ly = 5; ratiox = 3/4; ratioy = .5;
% Lx = 10; Ly = 5; ratiox = 1; ratioy = Ly/Lx;
[Nv, VX, VY, K, EToV] = unif_tri_mesh(round(ratiox*K1D),round(K1D*ratioy));
% VX = VX/max(abs(VX));  VY = VY/max(abs(VY));
% VX = (VX+1)*Lx; VY = VY*Ly;

% VX = VX + randn(size(VX));
StartUp2D;
BuildPeriodicMaps2D(max(VX)-min(VX),max(VY)-min(VY));

%% reference operators

% quadrature nodes
Nq = 3*N+2;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
[rq1D wq1D] = JacobiGQ(0,0,ceil(Nq/2));
% [rq1D wq1D] = JacobiGQ(0,0,N+1);
Vq = Vandermonde2D(N,rq,sq)/V;

% plotting nodes
[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
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
% DNr = diag(1./[wq;wfq])*QNr;
% DNs = diag(1./[wq;wfq])*QNs;

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

[r1 s1] = Nodes2D(Ngeo); [r1 s1] = xytors(r1,s1);
Vgeo = (Vandermonde2D(Ngeo,r,s)/Vandermonde2D(Ngeo,r1,s1))*(Vandermonde2D(N,r1,s1)/V);


% warp mesh near bdry
kk = 25;
w = .8;
Hx = (tanh(kk*(x+w))-tanh(kk*(x-w)))/2;
Hy = (tanh(kk*(y+w))-tanh(kk*(y-w)))/2;
H = 1-Hx.*Hy;
% H = 1; 

% warp mesh
kc = 3;
Lx = 1; Ly = 1;
x0 = 0; y0 = 0;
% dx = Lx*a*cos(kc*1/2*pi*(x-x0)/Lx).*cos(kc*3/2*pi*(y-y0)/Ly).*H;
% dy = Ly*a*sin(kc*3/2*pi*(x-x0)/Lx).*cos(kc*1/2*pi*(y-y0)/Ly).*H;
% x = x + dx;
% y = y + dy;

dy = Ly*a*cos(kc*3/2*pi*(x-x0)/Lx).*H;
dy = Vgeo*dy; % interp to Ngeo
y = y + dy;

dx = Lx*a*sin(kc*3/2*pi*(y-y0)/Ly).*H;
dx = Vgeo*dx; % interp to Ngeo
x = x + dx;

% x0 = 0; y0 = 0;
% x = x + Lx*a*cos(1/2*pi*(x-x0)/Lx).*cos(3/2*pi*(y-y0)/Ly);
% y = y + Ly*a*sin(pi*(x-x0)/Lx).*cos(1/2*pi*(y-y0)/Ly);

xq = Vq*x; yq = Vq*y;
xp = Vp*x; yp = Vp*y;
xf = Vf*x; yf = Vf*y;

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

if plotMesh
    rp1D = linspace(-1,1,100)';
    Vp1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,JacobiGL(0,0,N));
    Vfp = kron(eye(Nfaces),Vp1D);
    xfp = Vfp*x(Fmask(:),:);
    yfp = Vfp*y(Fmask(:),:);
    hold on
%     quiver(xf,yf,nxJ,nyJ)
    plot(xp,yp,'w.','markersize',16)
    plot(xfp,yfp,'k.')
    
%     L2err = nan;
axis off
axis equal
set(gca,'color','none')

% Ngeo = 1;
% [r1 s1] = Nodes2D(Ngeo); [r1 s1] = xytors(r1,s1);
% Vgeo = (Vandermonde2D(Ngeo,r,s)/Vandermonde2D(Ngeo,r1,s1))*(Vandermonde2D(N,r1,s1)/V);
% x = Vgeo*x;
% y = Vgeo*y;
% xfp = Vfp*x(Fmask(:),:);
% yfp = Vfp*y(Fmask(:),:);
% hold on
% plot(xfp,yfp,'r--','linewidth',2)
    return
end

%%

global avg
avg = @(x,y) .5*(x+y);

%% make phys frame ops

global wJq
wJq = diag(wq)*(J);

global Dx Dy
Dx = cell(K,1);
Dy = cell(K,1);
VqPN = cell(K,1);
VqLq = cell(K,1);
VpPq = cell(K,1);
E = cell(K,1);

if ~usePhysOps
    % FAKE operators 
    for e = 1:K
        [rxJ1, rxJ2] = meshgrid(rxJ(:,e));
        [sxJ1, sxJ2] = meshgrid(sxJ(:,e));
        [ryJ1, ryJ2] = meshgrid(ryJ(:,e));
        [syJ1, syJ2] = meshgrid(syJ(:,e));
        rxJK = avg(rxJ1,rxJ2);  sxJK = avg(sxJ1,sxJ2);
        ryJK = avg(ryJ1,ryJ2);  syJK = avg(syJ1,syJ2);
                
        Dx{e} = DNr.*rxJK + DNs.*sxJK;
        Dy{e} = DNr.*ryJK + DNs.*syJK;        
        MJ = Vq'*diag(wq.*J(:,e))*Vq;
        VqPN{e} = 2*Vq*(MJ\[Vq'*diag(wq) Vf'*diag(wfq)]);
        VqLq{e} = Vq*(MJ\(Vf'*diag(wfq)));        
        VpPq{e} = Vp*(MJ\(Vq'*diag(wq.*J(:,e))));        
        E{e} = Vf*(MJ\(Vq'*diag(wq.*J(:,e))));
    end
else
    
    % phys operators
    for e = 1:K
        
        Vphys = zeros(length(xq(:,e)),(N+1)*(N+2)/2);
        Vpphys = zeros(length(xp(:,e)),(N+1)*(N+2)/2);
        Vfphys = zeros(length(xf(:,e)),(N+1)*(N+2)/2);
        
        Vxphys = zeros(length(xq(:,e)),(N+1)*(N+2)/2);
        Vyphys = zeros(length(xq(:,e)),(N+1)*(N+2)/2);
        sk = 1;
        for i = 0:N
            for j = 0:N-i
                Vphys(:,sk) = xq(:,e).^i .* yq(:,e).^j;
                Vpphys(:,sk) = xp(:,e).^i .* yp(:,e).^j;
                Vfphys(:,sk) = xf(:,e).^i .* yf(:,e).^j;
                if i > 0
                    Vxphys(:,sk) = i*xq(:,e).^(i-1).*yq(:,e).^j;
                end
                if j > 0
                    Vyphys(:,sk) = xq(:,e).^i.*j.*yq(:,e).^(j-1);
                end
                sk = sk + 1;
            end
        end
        
        % Q: modes->qpts, Vphys: monomials->qpts
        % R = monomials->modes, inv(R): modes->monomials
        [Q R] = qr(Vphys,0);  % Q*R = Vphys
        Vq = Q;
        Vf = Vfphys/R;
        MJ = Vq'*diag(wq.*J(:,e))*Vq;
        Pq = MJ\(Vq'*diag(wq.*J(:,e))); % Pq: qpts -> modes
        
        Ee = Vf*Pq;
        E{e} = Ee;
        
        Dxmodal = Pq*(Vxphys/R);
        Dymodal = Pq*(Vyphys/R);
        Qx = diag(wq.*J(:,e))*Vq*Dxmodal*Pq;
        Qy = diag(wq.*J(:,e))*Vq*Dymodal*Pq;
        
        Bx = diag(nxJ(:,e).*wfq);
        By = diag(nyJ(:,e).*wfq);
        QNx = [Qx-.5*Ee'*Bx*Ee  .5*Ee'*Bx;
            -.5*Bx*Ee .5*Bx];
        QNy = [Qy-.5*Ee'*By*Ee .5*Ee'*By;
            -.5*By*Ee .5*By];
                
        Dx{e} = diag(1./[wq; wfq])*(QNx-QNx')*.5;
        Dy{e} = diag(1./[wq; wfq])*(QNy-QNy')*.5;
%         Dx{e} = diag(1./[wq; wfq])*(QNx);
%         Dy{e} = diag(1./[wq; wfq])*(QNy);
        VqPN{e} = 2*Vq*(MJ\[Vq'*diag(wq) Vf'*diag(wfq)]);
        VqLq{e} = Vq*(MJ\(Vf'*diag(wfq)));
        VpPq{e} = (Vpphys/R)*Pq;
        
%         % testing
%         u = xq(:,e).^N + yq(:,e).^N;
%         uf = Ee*u;
%         u'*diag(wJq(:,e))*VqPN{e}*Dx{e}*[u;uf]
        
    end
    
        
end


%%

uex = @(x,y) exp(-10*(x.^2+y.^2));
% uex = @(x,y) exp(-10*(x.^2+y.^2)).*sin(pi*x).*sin(pi*y);
x0 = 0;
u = VqPq*uex(xq-x0,yq);

initialErr = sqrt(sum(sum(wJq.*(uex(xq-x0,yq) - u).^2)))

% Runge-Kutta residual storage
res = zeros(Nq,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
CNh = max(CN*max(sJ(:)./Jf(:)));
dt = CFL*2/CNh;

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

figure(1)
for i = 1:Nsteps
    for INTRK = 1:5
        
        uM = zeros(size(nxJ));
        for e = 1:K
            uM(:,e) = E{e}*u(:,e);
        end
        
        % extra LF flux info
        rhs  = RHS2Dsimple(u,uM); 
        
        if INTRK==5
            rhstest(i) = sum(sum(u.*wJq.*rhs));
        end
        
        res = rk4a(INTRK)*res + dt*rhs;        
        u  = u  + rk4b(INTRK)*res;        
                
    end;
    
    if  (mod(i,5)==0 || i==Nsteps)
        clf
        vv = zeros(size(xp));
        for e = 1:K
            vv(:,e) = VpPq{e}*u(:,e);
        end
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f, N = %d, K1D = %d, rhstest = %g',dt*i,N,K1D,rhstest(i)))
        %         view(3)
        drawnow
    end
end

err = sqrt(sum(sum(wJq.*(u - uex(xq,yq)).^2)))

function [rhs] = RHS2Dsimple(uq,uM)

Globals2D;

global rxJ sxJ ryJ syJ nxJ nyJ
global mapPq
global fxS1 fyS1 fxS2 fyS2 fxS3 fyS3 avg
global DNr DNs VqPq VqLq VqPN

global Dx Dy

uP = uM(mapPq);

global tau

fSf = nxJ.*(uP) - tau*(uP-uM).*sJ;

uN = [uq; uM];
rhs = zeros(size(uq));
for e = 1:K
        
    % operators (replace w/phys frame)
    Dxe = Dx{e};
%     Dye = Dy{e};
    VqPNe = VqPN{e};
    VqLqe = VqLq{e};    
    
    divF = Dxe*uN(:,e);    
    rhs(:,e) = -(VqPNe*divF + VqLqe*fSf(:,e));    
end


end
