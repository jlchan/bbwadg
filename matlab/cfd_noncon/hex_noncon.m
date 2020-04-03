%% set up reference hex
clear
N = 1;
CFL = .25;
FinalTime = .33;

% interpolation nodes for geometric factors
[r1D w1D] = JacobiGL(0,0,N);
[r s t] = meshgrid(r1D);
r = r(:); s = s(:); t = t(:);
[wr ws wt] = meshgrid(w1D);
w = wr(:).*ws(:).*wt(:);

% quadrature nodes for collocation
[rq1D wq1D] = JacobiGQ(0,0,N);
% rq1D = rq1D*(1-1e-10);
[rq sq tq] = meshgrid(rq1D);
rq = rq(:); sq = sq(:); tq = tq(:);
[wr ws wt] = meshgrid(wq1D);
wq = wr(:).*ws(:).*wt(:);

% 1D matrices
V1D = Vandermonde1D(N,r1D);
D1D = GradVandermonde1D(N,r1D)/V1D;
I = eye(length(r1D));
Vq1D = Vandermonde1D(N,rq1D)/V1D;
M1D = Vq1D'*diag(wq1D)*Vq1D;
Pq1D = M1D\(Vq1D'*diag(wq1D));
Vf1D = Vandermonde1D(N,[-1;1])/Vandermonde1D(N,rq1D); % interp from quad to face points

% hybridized SBP matrices
Q1D = Pq1D'*M1D*D1D*Pq1D;
E1D = Vf1D;
B1D = diag([-1,1]);
QN = [Q1D-Q1D' E1D'*B1D;
    -B1D*E1D 0*B1D];
VN = [eye(N+1);Vf1D];
invWVNT = diag(1./wq1D)*VN';

% 3D matrices
V = kron(kron(V1D,V1D),V1D);
Vq = kron(kron(Vq1D,Vq1D),Vq1D);
M = Vq'*diag(wq)*Vq;
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));
Dr = kron(kron(I,D1D),I);
Ds = kron(kron(I,I),D1D);
Dt = kron(kron(D1D,I),I);

% points on a quad
[r2 s2] = meshgrid(rq1D);
[wr2 ws2] = meshgrid(wq1D);
r2 = r2(:); s2 = s2(:);
e = ones(size(r2));
wfq = wr2(:).*ws2(:);

% set face points according to r,s,t
rf = [-e; e; r2; r2; r2; r2];
sf = [r2; r2; -e; e; s2; s2];
tf = [s2; s2; s2; s2; -e; e];
wfq = repmat(wfq,6,1);

Vf = Vandermonde3DHex(N,rf,sf,tf)/Vandermonde3DHex(N,r,s,t); % interp nodes -> face nodes

% reference normals
nrJ = zeros((N+1)^2,6);
nsJ = zeros((N+1)^2,6);
ntJ = zeros((N+1)^2,6);
nrJ(:,1) = -1; nrJ(:,2) = 1;
nsJ(:,3) = -1; nsJ(:,4) = 1;
ntJ(:,5) = -1; ntJ(:,6) = 1;
nrJ = nrJ(:);
nsJ = nsJ(:);
ntJ = ntJ(:);

%% compute face indices for face points on "lines" of nodes

fidr = zeros(2,(N+1)^2);
fids = zeros(2,(N+1)^2);
fidt = zeros(2,(N+1)^2);
idr = (1:(N+1):(N+1)^2)';
ids = (1:N+1)';
idt = (1:(N+1)^2:(N+1)^3)';
for ii = 1:(N+1)^2
    rfr = Vf1D*rq(idr);    sfr = Vf1D*sq(idr);    tfr = Vf1D*tq(idr);
    rfs = Vf1D*rq(ids);    sfs = Vf1D*sq(ids);    tfs = Vf1D*tq(ids);
    rft = Vf1D*rq(idt);    sft = Vf1D*sq(idt);    tft = Vf1D*tq(idt);
    
    [rf1 rf2] = meshgrid(rf,rfr);
    [sf1 sf2] = meshgrid(sf,sfr);
    [tf1 tf2] = meshgrid(tf,tfr);
    D = abs(rf1 - rf2) + abs(sf1 - sf2) + abs(tf1 - tf2);
    fidr(1,ii) = find(D(1,:) < 1e-10);
    fidr(2,ii) = find(D(2,:) < 1e-10);
    
    [rf1 rf2] = meshgrid(rf,rfs);
    [sf1 sf2] = meshgrid(sf,sfs);
    [tf1 tf2] = meshgrid(tf,tfs);
    D = abs(rf1 - rf2) + abs(sf1 - sf2) + abs(tf1 - tf2);
    fids(1,ii) = find(D(1,:) < 1e-10);
    fids(2,ii) = find(D(2,:) < 1e-10);
    
    [rf1 rf2] = meshgrid(rf,rft);
    [sf1 sf2] = meshgrid(sf,sft);
    [tf1 tf2] = meshgrid(tf,tft);
    D = abs(rf1 - rf2) + abs(sf1 - sf2) + abs(tf1 - tf2);
    fidt(1,ii) = find(D(1,:) < 1e-10);
    fidt(2,ii) = find(D(2,:) < 1e-10);
    
    %     clf
    %     plot3(rq,sq,tq,'o')
    %     hold on
    %     plot3(rfr,sfr,tfr,'bx')
    %     plot3(rfs,sfs,tfs,'rx')
    %     plot3(rft,sft,tft,'kx')
    %     plot3(rq(idr),sq(idr),tq(idr),'bx-','linewidth',2) % moves in s
    %     plot3(rq(ids),sq(ids),tq(ids),'rs-','linewidth',2) % moves in r
    %     plot3(rq(idt),sq(idt),tq(idt),'m^-','linewidth',2) % moves in s
    %     pause
    
    % move onto next "line" of nodes
    idr = idr + 1;
    if mod(ii,N+1)==0
        idr = idr + (N+1)^2-(N+1);
    end
    ids = ids + (N+1);
    idt = idt + 1;
end

%% hybridized SBP - most operators for testing only

Qr = Pq'*M*Dr*Pq;
Qs = Pq'*M*Ds*Pq;
Qt = Pq'*M*Dt*Pq;

Ef = Vf*Pq;

Br = diag(wfq.*nrJ);
Bs = diag(wfq.*nsJ);
Bt = diag(wfq.*ntJ);

QNr = [Qr-Qr' Ef'*Br;
    -Br*Ef 0*Br];
QNs = [Qs-Qs' Ef'*Bs;
    -Bs*Ef 0*Bs];
QNt = [Qt-Qt' Ef'*Bt;
    -Bt*Ef 0*Bt];

% used for surface flux computations
Vfq = Vandermonde3DHex(N,rf,sf,tf)/Vandermonde3DHex(N,rq,sq,tq); % quadrature nodes -> face nodes
Vfq(abs(Vfq)<1e-8) = 0;
invWVfTW = diag(1./wq)*sparse(Vfq)'*diag(wfq);

VN = [eye((N+1)^3); Vfq];
invWVNTr = diag(1./wq)*VN';

% norm(QNr+QNr' - 2*blkdiag(0*Qr,Br),'fro')
% norm(sum(QNr,2))
% norm(sum(QNs,2))
% norm(sum(QNt,2))
% return

%% make physical elements

K = 2; 
x = zeros((N+1)^3,K);
y = zeros((N+1)^3,K);
z = zeros((N+1)^3,K);

% make multi-block "mesh"
x(:,1) = .5*(1+r)-1;
y(:,1) = s;
z(:,1) = t;
x(:,2) = .5*(1+r);
y(:,2) = s;
z(:,2) = t;

% form face node connectivities
xf = Vf*x;
yf = Vf*y;
zf = Vf*z;

Nfp = (N+1)^2;
Nfaces = 6;
mapM1 = reshape(1:Nfp*Nfaces,Nfp,Nfaces);
mapM2 = reshape(1:Nfp*Nfaces,Nfp,Nfaces) + Nfp*Nfaces;
mapP1 = mapM1;
mapP2 = mapM2;

% element 1
mapP1(:,1) = mapM2(:,2);
mapP1(:,2) = mapM2(:,1);
mapP1(:,3) = mapM1(:,4);
mapP1(:,4) = mapM1(:,3);
mapP1(:,5) = mapM1(:,6);
mapP1(:,6) = mapM1(:,5);

% element 2
mapP2(:,1) = mapM1(:,2);
mapP2(:,2) = mapM1(:,1);
mapP2(:,3) = mapM2(:,4);
mapP2(:,4) = mapM2(:,3);
mapP2(:,5) = mapM2(:,6);
mapP2(:,6) = mapM2(:,5);

% mapM = [mapM1 mapM2];
mapP = [mapP1 mapP2];
mapP = reshape(mapP,Nfp*Nfaces,K);

if 0
    for i = 1:length(mapP(:))
        plot3(xf,yf,zf,'o','linewidth',2,'markersize',16)
        hold on
        plot3(xf(i),yf(i),zf(i),'x','linewidth',2,'markersize',16)
        idP = mapP(i);
        plot3(xf(idP),yf(idP),zf(idP),'*','linewidth',2,'markersize',16)
        hold off
        view(-5,5)
        pause
    end
    text(xf(:)+.05,yf(:),zf(:),num2str((1:length(xf(:)))'))
end

%% compute geometric factors

% apply curved warping
a = 0*.125;
dx = cos(pi/2*x).*sin(pi*y).*sin(pi*z);
dy = sin(pi*x).*cos(pi/2*y).*sin(pi*z);
dz = sin(pi*x).*sin(pi*y).*cos(pi/2*z);
x = x + a*dx;
y = y + a*dy;
z = z + a*dz;

rxJ = Dt*((Ds*y).*z) - Ds*((Dt*y).*z);
sxJ = Dr*((Dt*y).*z) - Dt*((Dr*y).*z);
txJ = Ds*((Dr*y).*z) - Dr*((Ds*y).*z);

ryJ = -(Dt*((Ds*x).*z) - Ds*((Dt*x).*z));
syJ = -(Dr*((Dt*x).*z) - Dt*((Dr*x).*z));
tyJ = -(Ds*((Dr*x).*z) - Dr*((Ds*x).*z));

rzJ = -(Dt*((Ds*y).*x) - Ds*((Dt*y).*x));
szJ = -(Dr*((Dt*y).*x) - Dt*((Dr*y).*x));
tzJ = -(Ds*((Dr*y).*x) - Dr*((Ds*y).*x));

% compute normals
nxJ = nrJ.*(Vf*rxJ) + nsJ.*(Vf*sxJ) + ntJ.*(Vf*txJ);
nyJ = nrJ.*(Vf*ryJ) + nsJ.*(Vf*syJ) + ntJ.*(Vf*tyJ);
nzJ = nrJ.*(Vf*rzJ) + nsJ.*(Vf*szJ) + ntJ.*(Vf*tzJ);
sJ = sqrt(nxJ.^2 + nyJ.^2 + nzJ.^2);

% compute J
xr = Dr*x; xs = Ds*x; xt = Dt*x;
yr = Dr*y; ys = Ds*y; yt = Dt*y;
zr = Dr*z; zs = Ds*z; zt = Dt*z;
J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);

fprintf('GCL for elem  = %g\n',norm(Dr*rxJ + Ds*sxJ + Dt*txJ,'fro'))
fprintf('GCL for elem  = %g\n',norm(Dr*ryJ + Ds*syJ + Dt*tyJ,'fro'))
fprintf('GCL for elem  = %g\n',norm(Dr*rzJ + Ds*szJ + Dt*tzJ,'fro'))

xq = Vq*x;
yq = Vq*y;
zq = Vq*z;


%%

global gamma avg pfun betafun
gamma = 1.4;
avg = @(uL,uR) .5*(uL+uR);
pfun = @(rho, u, v, w, E) ((gamma - 1.0) * (E - .5 * rho .* (u .* u + v .* v + w .* w)));
betafun = @(rho, u, v, w, E) (rho ./ (2.0 * pfun(rho, u, v, w, E)));

rhoe = @(rho,rhou,rhov,rhow,E) E - .5*(rhou.^2+rhov.^2+rhow.^2)./rho;
pcons = @(rho,rhou,rhov,rhow,E) (gamma-1)*rhoe(rho,rhou,rhov,rhow,E);
sU = @(rho,rhou,rhov,rhow,E) log((gamma-1)*rhoe(rho,rhou,rhov,rhow,E)./(rho.^gamma));

sV = @(V1,V2,V3,V4,V5) gamma - V1 + (V2.^2+V3.^2+V4.^2)./(2*V5);
rhoeV  = @(V1,V2,V3,V4,V5) ((gamma-1)./((-V5).^gamma)).^(1/(gamma-1)).*exp(-sV(V1,V2,V3,V4,V5)/(gamma-1));

V1 = @(rho,rhou,rhov,rhow,E) (-E + rhoe(rho,rhou,rhov,rhow,E).*(gamma + 1 - sU(rho,rhou,rhov,rhow,E)))./(rhoe(rho,rhou,rhov,rhow,E));
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

%% initial conditions

time = 0;
rho = 2 + sin(pi*xq - time); 
u = ones(size(xq));
v = zeros(size(xq));
w = zeros(size(xq));
p = ones(size(xq));

rhou = rho.*u;
rhov = rho.*v;
rhow = rho.*w;
E    = p/(gamma-1) + .5*rho.*(u.^2+v.^2);

%% code to run solver

dt = CFL/(N+1)^2;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

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
    2802321613138.0/2924317926251.0 ...
    1.0];

res1 = zeros((N+1)^3,K);
res2 = zeros((N+1)^3,K);
res3 = zeros((N+1)^3,K);
res4 = zeros((N+1)^3,K);
res5 = zeros((N+1)^3,K);

for i = 1:Nsteps
    for INTRK = 1:5
        
        u    = rhou./rho;
        v    = rhov./rho;
        w    = rhow./rho;
        beta = betafun(rho,u,v,w,E);
        
        % entropy variables projection
        q1f = Vfq*V1(rho,rhou,rhov,rhow,E);
        q2f = Vfq*V2(rho,rhou,rhov,rhow,E);
        q3f = Vfq*V3(rho,rhou,rhov,rhow,E);
        q4f = Vfq*V4(rho,rhou,rhov,rhow,E);
        q5f = Vfq*V5(rho,rhou,rhov,rhow,E);
        
        rhof  = U1(q1f,q2f,q3f,q4f,q5f);
        rhouf = U2(q1f,q2f,q3f,q4f,q5f);
        rhovf = U3(q1f,q2f,q3f,q4f,q5f);
        rhowf = U4(q1f,q2f,q3f,q4f,q5f);
        Ef    = U5(q1f,q2f,q3f,q4f,q5f);
        
        uf    = rhouf./rhof;
        vf    = rhovf./rhof;
        wf    = rhowf./rhof;
        betaf = betafun(rhof,uf,vf,wf,Ef);
        
        % compute numerical flux
        rhoP  = rhof(mapP);
        uP    = uf(mapP);
        vP    = vf(mapP);
        wP    = wf(mapP);
        betaP = betaf(mapP);
        
        [FxSf FySf FzSf] = fluxes(rhof,rhoP,uf,uP,vf,vP,wf,wP,betaf,betaP);
        
        rhs1 = invWVfTW*(nxJ.*FxSf{1} + nyJ.*FySf{1} + nzJ.*FzSf{1});
        rhs2 = invWVfTW*(nxJ.*FxSf{2} + nyJ.*FySf{2} + nzJ.*FzSf{2});
        rhs3 = invWVfTW*(nxJ.*FxSf{3} + nyJ.*FySf{3} + nzJ.*FzSf{3});
        rhs4 = invWVfTW*(nxJ.*FxSf{4} + nyJ.*FySf{4} + nzJ.*FzSf{4});
        rhs5 = invWVfTW*(nxJ.*FxSf{5} + nyJ.*FySf{5} + nzJ.*FzSf{5});
        
        % slow version (was used for testing)
        if 0
            for e = 1:K
                
                rho_l  = [rho(:,e);  rhof(:,e)];
                u_l    = [u(:,e);    uf(:,e)];
                v_l    = [v(:,e);    vf(:,e)];
                w_l    = [w(:,e);    wf(:,e)];
                beta_l = [beta(:,e); betaf(:,e)];
                
                [rhox rhoy] = meshgrid(rho_l);
                [ux uy] = meshgrid(u_l);
                [vx vy] = meshgrid(v_l);
                [wx wy] = meshgrid(w_l);
                [betax betay] = meshgrid(beta_l);
                [FxS FyS FzS] = fluxes(rhox,rhoy,ux,uy,vx,vy,wx,wy,betax,betay);
                
                % affine case
                QNx = QNr*rxJ(1,e) + QNs*sxJ(1,e) + QNt*txJ(1,e);
                QNy = QNr*ryJ(1,e) + QNs*syJ(1,e) + QNt*tyJ(1,e);
                QNz = QNr*rzJ(1,e) + QNs*szJ(1,e) + QNt*tzJ(1,e);
                
                r1 = sum(QNx.*FxS{1} + QNy.*FyS{1} + QNz.*FzS{1},2);
                r2 = sum(QNx.*FxS{2} + QNy.*FyS{2} + QNz.*FzS{2},2);
                r3 = sum(QNx.*FxS{3} + QNy.*FyS{3} + QNz.*FzS{3},2);
                r4 = sum(QNx.*FxS{4} + QNy.*FyS{4} + QNz.*FzS{4},2);
                r5 = sum(QNx.*FxS{5} + QNy.*FyS{5} + QNz.*FzS{5},2);
                
                rhs1(:,e) = rhs1(:,e) + invWVNTr*r1;
                rhs2(:,e) = rhs2(:,e) + invWVNTr*r2;
                rhs3(:,e) = rhs3(:,e) + invWVNTr*r3;
                rhs4(:,e) = rhs4(:,e) + invWVNTr*r4;
                rhs5(:,e) = rhs5(:,e) + invWVNTr*r5;
            end
        end
        
        % compute volume terms
        for e = 1:K
            
            % indices for lines of nodes
            idrv = (1:(N+1):(N+1)^2)';
            idsv = (1:N+1)';
            idtv = (1:(N+1)^2:(N+1)^3)';
            
            for ii = 1:(N+1)^2
                
                idrf = fidr(:,ii);
                idsf = fids(:,ii);
                idtf = fidt(:,ii);
                
                % ========= differentiate along r-lines
                
                rho_l  = [rho(idrv,e);  rhof(idrf,e)];
                u_l    = [u(idrv,e);    uf(idrf,e)];
                v_l    = [v(idrv,e);    vf(idrf,e)];
                w_l    = [w(idrv,e);    wf(idrf,e)];
                beta_l = [beta(idrv,e); betaf(idrf,e)];
                
                [rhox rhoy] = meshgrid(rho_l);
                [ux uy] = meshgrid(u_l);
                [vx vy] = meshgrid(v_l);
                [wx wy] = meshgrid(w_l);
                [betax betay] = meshgrid(beta_l);
                [FxS FyS FzS] = fluxes(rhox,rhoy,ux,uy,vx,vy,wx,wy,betax,betay);
                
                % affine case
                Fr1 = FxS{1}*rxJ(1,e) + FyS{1}*ryJ(1,e) + FzS{1}*rzJ(1,e);
                Fr2 = FxS{2}*rxJ(1,e) + FyS{2}*ryJ(1,e) + FzS{2}*rzJ(1,e);
                Fr3 = FxS{3}*rxJ(1,e) + FyS{3}*ryJ(1,e) + FzS{3}*rzJ(1,e);
                Fr4 = FxS{4}*rxJ(1,e) + FyS{4}*ryJ(1,e) + FzS{4}*rzJ(1,e);
                Fr5 = FxS{5}*rxJ(1,e) + FyS{5}*ryJ(1,e) + FzS{5}*rzJ(1,e);
                
                rhs1(idrv,e) = rhs1(idrv,e) + invWVNT*sum(QN.*Fr1,2);
                rhs2(idrv,e) = rhs2(idrv,e) + invWVNT*sum(QN.*Fr2,2);
                rhs3(idrv,e) = rhs3(idrv,e) + invWVNT*sum(QN.*Fr3,2);
                rhs4(idrv,e) = rhs4(idrv,e) + invWVNT*sum(QN.*Fr4,2);
                rhs5(idrv,e) = rhs5(idrv,e) + invWVNT*sum(QN.*Fr5,2);
                
                % ========= differentiate along s-lines
                
                rho_l  = [rho(idsv,e);  rhof(idsf,e)];
                u_l    = [u(idsv,e);    uf(idsf,e)];
                v_l    = [v(idsv,e);    vf(idsf,e)];
                w_l    = [w(idsv,e);    wf(idsf,e)];
                beta_l = [beta(idsv,e); betaf(idsf,e)];
                
                [rhox rhoy] = meshgrid(rho_l);
                [ux uy] = meshgrid(u_l);
                [vx vy] = meshgrid(v_l);
                [wx wy] = meshgrid(w_l);
                [betax betay] = meshgrid(beta_l);
                [FxS FyS FzS] = fluxes(rhox,rhoy,ux,uy,vx,vy,wx,wy,betax,betay);
                
                % affine case
                Fs1 = FxS{1}*sxJ(1,e) + FyS{1}*syJ(1,e) + FzS{1}*szJ(1,e);
                Fs2 = FxS{2}*sxJ(1,e) + FyS{2}*syJ(1,e) + FzS{2}*szJ(1,e);
                Fs3 = FxS{3}*sxJ(1,e) + FyS{3}*syJ(1,e) + FzS{3}*szJ(1,e);
                Fs4 = FxS{4}*sxJ(1,e) + FyS{4}*syJ(1,e) + FzS{4}*szJ(1,e);
                Fs5 = FxS{5}*sxJ(1,e) + FyS{5}*syJ(1,e) + FzS{5}*szJ(1,e);
                rhs1(idsv,e) = rhs1(idsv,e) + invWVNT*sum(QN.*Fs1,2);
                rhs2(idsv,e) = rhs2(idsv,e) + invWVNT*sum(QN.*Fs2,2);
                rhs3(idsv,e) = rhs3(idsv,e) + invWVNT*sum(QN.*Fs3,2);
                rhs4(idsv,e) = rhs4(idsv,e) + invWVNT*sum(QN.*Fs4,2);
                rhs5(idsv,e) = rhs5(idsv,e) + invWVNT*sum(QN.*Fs5,2);
                
                % ========= differentiate along t-lines
                
                rho_l  = [rho(idtv,e);  rhof(idtf,e)];
                u_l    = [u(idtv,e);    uf(idtf,e)];
                v_l    = [v(idtv,e);    vf(idtf,e)];
                w_l    = [w(idtv,e);    wf(idtf,e)];
                beta_l = [beta(idtv,e); betaf(idtf,e)];
                
                [rhox rhoy] = meshgrid(rho_l);
                [ux uy] = meshgrid(u_l);
                [vx vy] = meshgrid(v_l);
                [wx wy] = meshgrid(w_l);
                [betax betay] = meshgrid(beta_l);
                [FxS FyS FzS] = fluxes(rhox,rhoy,ux,uy,vx,vy,wx,wy,betax,betay);
                
                % affine case
                Ft1 = FxS{1}*txJ(1,e) + FyS{1}*tyJ(1,e) + FzS{1}*tzJ(1,e);
                Ft2 = FxS{2}*txJ(1,e) + FyS{2}*tyJ(1,e) + FzS{2}*tzJ(1,e);
                Ft3 = FxS{3}*txJ(1,e) + FyS{3}*tyJ(1,e) + FzS{3}*tzJ(1,e);
                Ft4 = FxS{4}*txJ(1,e) + FyS{4}*tyJ(1,e) + FzS{4}*tzJ(1,e);
                Ft5 = FxS{5}*txJ(1,e) + FyS{5}*tyJ(1,e) + FzS{5}*tzJ(1,e);
                rhs1(idtv,e) = rhs1(idtv,e) + invWVNT*sum(QN.*Ft1,2);
                rhs2(idtv,e) = rhs2(idtv,e) + invWVNT*sum(QN.*Ft2,2);
                rhs3(idtv,e) = rhs3(idtv,e) + invWVNT*sum(QN.*Ft3,2);
                rhs4(idtv,e) = rhs4(idtv,e) + invWVNT*sum(QN.*Ft4,2);
                rhs5(idtv,e) = rhs5(idtv,e) + invWVNT*sum(QN.*Ft5,2);
                
                % move onto next "line" of nodes
                idrv = idrv + 1;
                if mod(ii,N+1)==0
                    idrv = idrv + (N+1)^2-(N+1);
                end
                idsv = idsv + (N+1);
                idtv = idtv + 1;
                
            end
        end
        
        
        if (INTRK==5)
            % rhstest
            r1 = diag(wq)*rhs1;
            r2 = diag(wq)*rhs2;
            r3 = diag(wq)*rhs3;
            r4 = diag(wq)*rhs4;
            r5 = diag(wq)*rhs5;
            v1 = V1(rho,rhou,rhov,rhow,E);
            v2 = V2(rho,rhou,rhov,rhow,E);
            v3 = V3(rho,rhou,rhov,rhow,E);
            v4 = V4(rho,rhou,rhov,rhow,E);
            v5 = V5(rho,rhou,rhov,rhow,E);
            rhstest(i) = sum(sum(v1.*r1 + v2.*r2 + v3.*r3 + v4.*r4 + v5.*r5));
        end
        
        rhs1 = -rhs1./J;
        rhs2 = -rhs2./J;
        rhs3 = -rhs3./J;
        rhs4 = -rhs4./J;
        rhs5 = -rhs5./J;
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        res4 = rk4a(INTRK)*res4 + dt*rhs4;
        res5 = rk4a(INTRK)*res5 + dt*rhs5;
        
        rho  = rho  + rk4b(INTRK)*res1;
        rhou = rhou + rk4b(INTRK)*res2;
        rhov = rhov + rk4b(INTRK)*res3;
        rhow = rhow + rk4b(INTRK)*res4;
        E    = E    + rk4b(INTRK)*res5;
    end
    
    if mod(i,25)==0 || i==Nsteps
        clf
        h = color_line3(xq,yq,zq,rho,'.');
        set(h,'markersize',64)
        axis tight
        view(-5,5)
        set(gca,'fontsize',15)
        title(sprintf('at tstep = %d / %d, time = %f, rhstest = %g\n',i,Nsteps,dt*i,rhstest(i)))
        drawnow limitrate
    end
end

rhoex = 2 + sin(pi*(xq - FinalTime));

% figure
% subplot(1,2,1)
% h = color_line3(xq,yq,zq,rho,'.');
% set(h,'markersize',32)
% axis tight
% view(-5,5)
% title(sprintf('at tstep = %d / %d, time = %f\n',i,Nsteps,dt*i))
% colorbar
% set(gca,'fontsize',15)
% subplot(1,2,2)
% h = color_line3(xq,yq,zq,rhoex,'.');
% set(h,'markersize',32)
% axis tight
% view(-5,5)
% title(sprintf('at tstep = %d / %d, time = %f\n',i,Nsteps,dt*i))
% colorbar
% set(gca,'fontsize',15)
max(max(abs(rho - rhoex)))

%%

function [FxS FyS FzS] = fluxes(rhoL,rhoR,uL,uR,vL,vR,wL,wR,betaL,betaR)

global gamma
global avg pfun
rholog = logmean(rhoL, rhoR);
rhoavg = avg(rhoL, rhoR);
uavg = avg(uL, uR);
vavg = avg(vL, vR);
wavg = avg(wL, wR);
% vnavg = 2.0 * (uavg.^2 + vavg.^2 + wavg.^2) - ...
%     (avg(uL.^2, uR.^2) + avg(vL.^2, vR.^2) + avg(wL.^2, wR.^2));
vnavg = uL.*uR + vL.*vR + wL.*wR; % this is way faster than the above...
pa = rhoavg ./ (2.0 * avg(betaL, betaR));
f4aux = rholog ./ (2.0 * (gamma - 1.0) * logmean(betaL, betaR)) + pa + .5 * rholog .* vnavg;

FxS{1} = rholog.*uavg;
FyS{1} = rholog.*vavg;
FzS{1} = rholog.*wavg;

FxS{2} = FxS{1}.*uavg + pa;
FyS{2} = FyS{1}.*uavg;
FzS{2} = FzS{1}.*uavg;

FxS{3} = FxS{1}.*vavg;
FyS{3} = FyS{1}.*vavg + pa;
FzS{3} = FzS{1}.*vavg;

FxS{4} = FxS{1}.*wavg;
FyS{4} = FyS{1}.*wavg;
FzS{4} = FzS{1}.*wavg + pa;

FxS{5} = f4aux.*uavg;
FyS{5} = f4aux.*vavg;
FzS{5} = f4aux.*wavg;

end



