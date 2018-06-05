clear -globals
for NB = [4]
    sk = 1;
    for Ksub = [4 8 16 32]
        Ksub
        [L2err_wadg(sk) pwadg MM] = WaveQuad_IGA(NB,Ksub,1,0,1);
        [L2err_proj(sk) pproj MM] = WaveQuad_IGA(NB,Ksub,1,0,0);
        dp = pwadg-pproj;
        hh(sk) = 2/Ksub;
        difference(sk) = dp'*MM*dp;
        sk = sk + 1;
    end
    
    loglog(hh,L2err_proj,'o--')
    hold on
    loglog(hh,L2err_wadg,'x--')
    loglog(hh,difference,'s--')
   
    print_pgf_coordinates(hh,L2err_proj)
    print_pgf_coordinates(hh,L2err_wadg)
    print_pgf_coordinates(hh,difference)

end

function [L2err p MM] = WaveQuad_IGA(NBin,Ksubin,K1Din,smoothKnotsin,useWADGin)

clear global
Globals2D;

global NB Ksub

if nargin==0
    NB = 4;
    Ksub = 32;
    K1D = 2;
    smoothKnots = 25;

else
    NB = NBin;
    Ksub = Ksubin;
    K1D = K1Din;
    smoothKnots = smoothKnotsin;

end

a = .125;
a = .28;

N = NB+Ksub-1;
dofs = (N+1)^2*K1D^2;

FinalTime = .5;

% Read in Mesh
[Nv, VX, VY, K, EToV] = QuadMesh2D(K1D);

% a = .1;
% ids = find(abs(abs(VX)-1) > 1e-8 & abs(abs(VY)-1) > 1e-8);
% VX(ids) = VX(ids) + a;
% VY(ids) = VY(ids) + a;

StartUpQuad2D;
% plot(x,y,'o');return
% BuildPeriodicMaps2D(2,2);
% plot(x(vmapM),y(vmapM),'o')
% hold on
% plot(x(vmapP),y(vmapP),'x')

%%
global rxJ sxJ ryJ syJ nxq nyq Jq sJq
% global Vq Vrq Vsq Pq Prq Psq
global Vfq Pfq

% non-affine mappings
x = x + a*cos(.5*3*pi*y).*cos(pi/2*x);
y = y + a*sin(.5*3*pi*x).*cos(pi/2*y);

% plot(x,y,'o');return

% isoparametric - interpolate geofacs
[rx,sx,ry,sy,J] = GeometricFactorsQuad2D(x,y,Dr,Ds);
[nx, ny, sJ] = NormalsQuad2D();


%% initialize TP operators

rp1D = linspace(-1,1,150);
[rp sp] = meshgrid(rp1D);
rp = rp(:); sp = sp(:);
Vp1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,JacobiGL(0,0,N));
for e = 1:K
    xp(:,e) = reshape(Vp1D*reshape(x(:,e),N+1,N+1)*Vp1D',length(rp),1);
    yp(:,e) = reshape(Vp1D*reshape(y(:,e),N+1,N+1)*Vp1D',length(rp),1);
end

[r1D] = JacobiGL(0,0,N);
[~, ~, ~, ~, rq1D, wq1D] = bsplineVDM(NB,Ksub,r1D,smoothKnots); % VDM for interp, mass, M\S

[rq sq] = meshgrid(rq1D); rq = rq(:); sq = sq(:);
[wrq wsq] = meshgrid(wq1D); wq = wrq(:).*wsq(:);
V1D = Vandermonde1D(N,r1D);
Vq1D = Vandermonde1D(N,rq1D)/V1D;

Nq = length(rq1D);
rxj = (rx.*J); sxj = (sx.*J);
ryj = (ry.*J); syj = (sy.*J);

xq = zeros(Nq^2,K);
yq = zeros(Nq^2,K);
Jq = zeros(Nq^2,K);
rxJ = zeros(Nq^2,K); sxJ = zeros(Nq^2,K);
ryJ = zeros(Nq^2,K); syJ = zeros(Nq^2,K);
for e = 1:K
    xq(:,e) = reshape(Vq1D*reshape(x(:,e),N+1,N+1)*Vq1D',Nq^2,1);
    yq(:,e) = reshape(Vq1D*reshape(y(:,e),N+1,N+1)*Vq1D',Nq^2,1);
    Jq(:,e) = reshape(Vq1D*reshape(J(:,e),N+1,N+1)*Vq1D',Nq^2,1);
    rxJ(:,e) = reshape(Vq1D*reshape(rxj(:,e),N+1,N+1)*Vq1D',Nq^2,1);
    sxJ(:,e) = reshape(Vq1D*reshape(sxj(:,e),N+1,N+1)*Vq1D',Nq^2,1);
    ryJ(:,e) = reshape(Vq1D*reshape(ryj(:,e),N+1,N+1)*Vq1D',Nq^2,1);
    syJ(:,e) = reshape(Vq1D*reshape(syj(:,e),N+1,N+1)*Vq1D',Nq^2,1);
end
wJq = spdiag(wq)*Jq;

Vfq = blkdiag(Vq1D,Vq1D,Vq1D,Vq1D);
nxq = Vfq*nx; nyq = Vfq*ny;
sJq = Vfq*sJ;

% return

global Vq1D Vrq1D Pq1D Prq1D 
% if Ksub==1
%     V1D = Vandermonde1D(N,r1D);
%     D1D = GradVandermonde1D(N,r1D)/V1D;
%     Vp1D = Vandermonde1D(N,rp1D)/V1D;
%     Vq1D = Vandermonde1D(N,rq1D)/V1D;    
%     Vrq1D = Vq1D*D1D;
%     M1D = Vq1D'*diag(wq1D)*Vq1D;
%     invM1D = V1D*V1D';
%     Pq1D = invM1D*Vq1D'*diag(wq1D);
% else
    [BVDM M1D D1D R, ~, ~, Vq1D Vrq1D] = bsplineVDM(NB,Ksub,r1D,smoothKnots); % VDM for interp, mass, M\S    
    Vp1D = bsplineVDM(NB,Ksub,rp1D,smoothKnots);
    
    invM1D = inv(M1D);
    Pq1D = invM1D*Vq1D'*diag(wq1D);
    Prq1D = invM1D*Vrq1D'*diag(wq1D);        
% end
% Vp = kron(Vp1D,Vp1D);

Mf = zeros(N+1,4*(N+1));
for f = 1:4
    Mf(Fmask(:,f),(1:N+1) + (f-1)*(N+1)) = M1D;
end
for i = 1:size(LIFT,2)
    LIFT(:,i) = reshape(invM1D * reshape(Mf(:,i),N+1,N+1)*invM1D',(N+1)^2,1);
end
LIFT(abs(LIFT)<1e-8) = 0;

% project onto face space then apply usual lift
Vfq = blkdiag(Vq1D,Vq1D,Vq1D,Vq1D);
Pfq = LIFT*kron(eye(Nfaces),invM1D)*Vfq'*diag(repmat(wq1D,Nfaces,1));

% global Grr Grs Gss Gst
% Grr = (rxJ.^2 + ryJ.^2)./(Jq.^2);
% Grs = (rxJ.*sxJ + ryJ.*syJ)./(Jq.^2);
% Gss = (sxJ.^2 + syJ.^2)./(Jq.^2);


%% check eigs

if 0 && nargin==0 && Np*K < 3000
    U = zeros(Np*K,1);
    A = zeros(Np*K);
    for i = 1:Np*K
        U(i) = 1;
        p = reshape(U(:,1),Np,K);
        [rhs] = acousticsRHS2Dq(p);
        A(:,i) = rhs(:);
        U(i) = 0;
        if (mod(i,ceil(3*Np*K/10))==0)
            disp(sprintf('on i = %d out of %d\n',i,3*Np*K))
        end
    end
    lam = eig(A);
    plot(lam,'o');
    title(sprintf('largest real part = %g\n',max(real(lam))))
    max(abs(lam))
    keyboard
end

%% global mass matrix

global useWADG MJ Mref
useWADG = useWADGin;
Vq = kron(sparse(Vq1D),sparse(Vq1D));
MJ = cell(K,1);
Mref = Vq'*spdiag(wq)*Vq;
for e = 1:K
    MJ{e} = Vq'*spdiag(wJq(:,e))*Vq;    
end
MM = blkdiag(MJ{:});

%% estimate timestep

CT = max((NB+1)*Ksub,(NB+1)^2);
h = min(min(Jq)./max(sJq));
dt = h / CT;

%% initial cond

% x0 = 0; y0 = .1;
% pex = @(x,y,t) exp(-15^2*((x-x0).^2 + (y-y0).^2));
k = 3;
pex = @(x,y,t) cos(k*pi*x/2).*cos(k*pi*y/2).*cos(sqrt(2)*.5*k*pi*t);

% p = Pq*pex(xq,yq,0);
% for e = 1:K
%     p(:,e) = (Vq'*diag(wq.*Jq(:,e))*Vq)\(Vq'* (wq.*Jq(:,e).*pex(xq(:,e),yq(:,e),0)));
% end
Nq = length(rq1D);
for e = 1:K
    p(:,e) = reshape(Pq1D*reshape(pex(xq(:,e),yq(:,e),0),Nq,Nq)*Pq1D',(N+1)^2,1);
end

u = zeros(Np, K);

pq = zeros(Nq^2,K);
for e = 1:K
    pK = Vq1D*reshape(p(:,e),N+1,N+1)*Vq1D';
    pq(:,e) = pK(:);
end
if nargin==0
    err = wJq.*(pq - pex(xq,yq,0)).^2;
    init_cond_err = sqrt(sum(err(:)))
end

% keyboard
%%

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K); resp = zeros(Np,K);

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

% outer time step loop
for tstep = 1:Nsteps
    for INTRK = 1:5
        
        timelocal = tstep*dt + rk4c(INTRK)*dt;
        [rhsu] = acousticsRHS2Dq(p);
        rhsp = u;
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        p = p+rk4b(INTRK)*resp;
        
    end;
    
    if 1 && nargin==0 && (mod(tstep,5)==0 || tstep==Nsteps)
        clf
        vv = zeros(length(rp),K);
        for e = 1:K
            vv(:,e) = reshape(Vp1D*reshape(p(:,e),N+1,N+1)*Vp1D',length(rp),1);
        end
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f',tstep*dt))
        drawnow
    end
    
    if nargin==0 && mod(tstep,25)==0
        disp(sprintf('on tstep %d out of %d\n',tstep,Nsteps))
    end
end


pq = zeros(Nq^2,K);
for e = 1:K
    pK = Vq1D*reshape(p(:,e),N+1,N+1)*Vq1D';
    pq(:,e) = pK(:);
end

err = wJq.*(pq - pex(xq,yq,FinalTime)).^2;
L2err = sqrt(sum(err(:)));
end

function [rhs] = acousticsRHS2Dq(p)

Globals2D;

global rxJ sxJ ryJ syJ nxq nyq Jq sJq
% global Vq Vrq Vsq Pq Prq Psq
global Vfq Pfq
global Vq1D Vrq1D Pq1D Prq1D 

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapM) - p(vmapP);
dp(mapB) = 2*p(vmapB);
dpq = Vfq*dp;

% compute q
% prq = Vrq*p;
% psq = Vsq*p;
Nq = size(Vq1D,1);
prq = zeros(Nq^2,K);
psq = zeros(Nq^2,K);
for e = 1:K
    prq(:,e) = reshape(Vrq1D*reshape(p(:,e),N+1,N+1)*Vq1D',Nq^2,1);
    psq(:,e) = reshape(Vq1D*reshape(p(:,e),N+1,N+1)*Vrq1D',Nq^2,1);
end

%ux = Pq*((Vq*Pq*(rxJ.*prq + sxJ.*psq))./Jq);
%uy = Pq*((Vq*Pq*(ryJ.*prq + syJ.*psq))./Jq);
pxq = (rxJ.*prq + sxJ.*psq)./Jq;
pyq = (ryJ.*prq + syJ.*psq)./Jq;
px = zeros(Np,K);
py = zeros(Np,K);
for e = 1:K    
    px(:,e) = reshape(Pq1D*reshape(pxq(:,e),Nq,Nq)*Pq1D',Np,1);
    py(:,e) = reshape(Pq1D*reshape(pyq(:,e),Nq,Nq)*Pq1D',Np,1);
end

qx = px - .5*Pfq*(dpq.*nxq.*sJq);
qy = py - .5*Pfq*(dpq.*nyq.*sJq);

qxr = zeros(Nq^2,K); qxs = zeros(Nq^2,K);
qyr = zeros(Nq^2,K); qys = zeros(Nq^2,K);
for e = 1:K
    qxr(:,e) = reshape(Vrq1D*reshape(qx(:,e),N+1,N+1)*Vq1D',Nq^2,1);
    qxs(:,e) = reshape(Vq1D*reshape(qx(:,e),N+1,N+1)*Vrq1D',Nq^2,1);
    qyr(:,e) = reshape(Vrq1D*reshape(qy(:,e),N+1,N+1)*Vq1D',Nq^2,1);
    qys(:,e) = reshape(Vq1D*reshape(qy(:,e),N+1,N+1)*Vrq1D',Nq^2,1);
end
% qxr = Vrq*qx; qxs = Vsq*qx;
% qyr = Vrq*qy; qys = Vsq*qy;
divQ = rxJ.*qxr + sxJ.*qxs + ryJ.*qyr + syJ.*qys;

dqx = Vfq*reshape(qx(vmapM)-.5*(px(vmapM)+px(vmapP)),Nfp*Nfaces,K);
dqy = Vfq*reshape(qy(vmapM)-.5*(py(vmapM)+py(vmapP)),Nfp*Nfaces,K);

global NB Ksub
CT = max((NB+1)*Ksub,(NB+1)^2);
hmin = min(min(Jq)./max(sJq));
tau = CT / hmin;

dqn = nxq.*dqx + nyq.*dqy;
fluxq = dqn + tau*dpq;

rhsflux = Pfq*(sJq.*fluxq);
% rhs = Pq*((Vq*rhs)./Jq);

global useWADG MJ Mref
rhs = zeros(Np,K);
for e = 1:K
    rp = Pq1D*reshape(divQ(:,e),Nq,Nq)*Pq1D';
    rhsp = rp(:) - rhsflux(:,e);
    if useWADG
        JqK = reshape(Jq(:,e),Nq,Nq);
        rhs(:,e) = reshape(Pq1D*((Vq1D*reshape(rhsp,N+1,N+1)*Vq1D')./JqK)*Pq1D',Np,1);
    else
        rhs(:,e) = MJ{e}\(Mref*rhsp);
    end
end
end

