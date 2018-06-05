clear -globals
clear
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


function [L2err p MM] = WaveQuad_IGA(NB,Ksub,K1D,smoothKnots,useWADGin)


Globals2D;
if nargin==0
    NB = 4;
    Ksub = 32;
    K1D = 1;
    smoothKnots = 75;
end

a = .125;
a = .28;

N = NB+Ksub-1;
dofs = (N+1)^2*K1D^2;

FinalTime = .1;

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
global Pfq Vfq

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
% Vp = Vandermonde2DQuad(N,rp,sp)/V;
Vp1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,JacobiGL(0,0,N));
for e = 1:K
    xp(:,e) = reshape(Vp1D*reshape(x(:,e),N+1,N+1)*Vp1D',length(rp),1);
    yp(:,e) = reshape(Vp1D*reshape(y(:,e),N+1,N+1)*Vp1D',length(rp),1);
end

[r1D] = JacobiGL(0,0,N);
[~, ~, ~, ~, rq1D, wq1D] = bsplineVDM(NB,Ksub,r1D,smoothKnots); % VDM for interp, mass, M\S



global Vq1D Vrq1D Pq1D Prq1D
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

Vfq = sparse(blkdiag(Vq1D,Vq1D,Vq1D,Vq1D));
nxq = Vfq*nx; nyq = Vfq*ny;
sJq = Vfq*sJ;


if Ksub==1
    V1D = Vandermonde1D(N,r1D);
    D1D = GradVandermonde1D(N,r1D)/V1D;
    Vp1D = Vandermonde1D(N,rp1D)/V1D;
    Vq1D = Vandermonde1D(N,rq1D)/V1D;
    M1D = Vq1D'*diag(wq1D)*Vq1D;
    invM1D = V1D*V1D';
else
    [BVDM M1D D1D R, ~, ~, Vq1D Vrq1D] = bsplineVDM(NB,Ksub,r1D,smoothKnots); % VDM for interp, mass, M\S
    Vp1D = bsplineVDM(NB,Ksub,rp1D,smoothKnots);
    
    invM1D = inv(M1D);
    Pq1D = invM1D*Vq1D'*diag(wq1D);
    Prq1D = invM1D*Vrq1D'*diag(wq1D);
    
end

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

% keyboard

%% check eigs

if 0 && nargin==0 && 3*Np*K < 3000
    U = zeros(Np*K,3);
    A = zeros(3*Np*K);
    for i = 1:3*Np*K
        U(i) = 1;
        p = reshape(U(:,1),Np,K);
        u = reshape(U(:,2),Np,K);
        v = reshape(U(:,3),Np,K);
        [rhsp, rhsu, rhsv] = acousticsRHS2Dq(p,u,v);
        A(:,i) = [rhsp(:);rhsu(:);rhsv(:)];
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


%%

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

U = randn((N+1)^2*K,3);
for i = 1:10
    Uprev = U;
    [rhsp, rhsu, rhsv] = acousticsRHS2Dq(reshape(U(:,1),Np,K),reshape(U(:,2),Np,K),reshape(U(:,3),Np,K));
    U(:,1) = rhsp(:);
    U(:,2) = rhsu(:);
    U(:,3) = rhsv(:);
    
    lam = Uprev(:)'*U(:) / norm(Uprev(:));
    U = U/norm(U(:));
end
dt = .5/abs(lam);



%% initial cond

% x0 = 0; y0 = .1;
% pex = @(x,y,t) exp(-15^2*((x-x0).^2 + (y-y0).^2));
k = 3;
pex = @(x,y,t) cos(k*pi*x/2).*cos(k*pi*y/2).*cos(sqrt(2)*.5*k*pi*t);

for e = 1:K
    p(:,e) = reshape(Pq1D*reshape(pex(xq(:,e),yq(:,e),0),Nq,Nq)*Pq1D',(N+1)^2,1);
end

% p = VDM\pex(x,y,0);
u = zeros(Np, K);
v = zeros(Np, K);

pq = zeros(Nq^2,K);
for e = 1:K
    pK = Vq1D*reshape(p(:,e),N+1,N+1)*Vq1D';
    pq(:,e) = pK(:);
end
if nargin==0
    err = wJq.*(pq - pex(xq,yq,0)).^2;
    init_cond_err = sqrt(sum(err(:)))
end

% return
% keyboard
%%

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

% outer time step loop
for tstep = 1:Nsteps
    for INTRK = 1:5
        
        timelocal = tstep*dt + rk4c(INTRK)*dt;
        [rhsp, rhsu, rhsv] = acousticsRHS2Dq(p,u,v);
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
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

function [rhsp, rhsu, rhsv] = acousticsRHS2Dq(p,u,v)

Globals2D;

global rxJ sxJ ryJ syJ nxq nyq Jq sJq
% global Vq Vrq Vsq Pq Prq Psq
global Vfq Pfq
global Vq1D Vrq1D Pq1D Prq1D

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

uavg = zeros(Nfp*Nfaces,K); uavg(:) = u(vmapP)+u(vmapM);
vavg = zeros(Nfp*Nfaces,K); vavg(:) = v(vmapP)+v(vmapM);

% uavg(mapB) = u(vmapB);
% vavg(mapB) = v(vmapB);
% du(mapB) = 0;
% dv(mapB) = 0;
dp(mapB) = -2*p(vmapB);

dp = Vfq*dp;
du = Vfq*du;
dv = Vfq*dv;
uavg = Vfq*uavg;
vavg = Vfq*vavg;

% Impose reflective boundary conditions (p+ = -p-)
ndotUavg = nxq.*uavg + nyq.*vavg;
ndotdU = nxq.*du + nyq.*dv;

tau = 1;
fluxp =  tau*dp - ndotUavg;
fluxu =  (tau*ndotdU - dp).*nxq;
fluxv =  (tau*ndotdU - dp).*nyq;

% pr = Vrq*p;
% ps = Vsq*p;
Nq = size(Vrq1D,1);
pr = zeros(Nq^2,K); ps = zeros(Nq^2,K);
uq = zeros(Nq^2,K); vq = zeros(Nq^2,K);
for e = 1:K
    prK = Vrq1D*reshape(p(:,e),N+1,N+1)*Vq1D';
    psK = Vq1D*reshape(p(:,e),N+1,N+1)*Vrq1D';
    
    uqK = Vq1D*reshape(u(:,e),N+1,N+1)*Vq1D';
    vqK = Vq1D*reshape(v(:,e),N+1,N+1)*Vq1D';
    
    uq(:,e) = uqK(:);
    vq(:,e) = vqK(:);
    pr(:,e) = prK(:);
    ps(:,e) = psK(:);
end
px = rxJ.*pr + sxJ.*ps;
py = ryJ.*pr + syJ.*ps;

% dpdx = Pq*(px);
% dpdy = Pq*(py);
Ur = rxJ.*uq + ryJ.*vq;
Us = sxJ.*uq + syJ.*vq;
% divU = -(Prq*(Ur) + Psq*(Us));
for e = 1:K
    pxK = Pq1D*reshape(px(:,e),Nq,Nq)*Pq1D';
    pyK = Pq1D*reshape(py(:,e),Nq,Nq)*Pq1D';
    divUK = Prq1D*reshape(Ur(:,e),Nq,Nq)*Pq1D' + Pq1D*reshape(Us(:,e),Nq,Nq)*Prq1D';
    dpdx(:,e) = pxK(:);
    dpdy(:,e) = pyK(:);
    divU(:,e) = -divUK(:);
end

% compute right hand sides of the PDE's
rhsp =  (-divU) + Pfq*(sJq.*fluxp)/2.0;
rhsu =  (-dpdx) + Pfq*(sJq.*fluxu)/2.0;
rhsv =  (-dpdy) + Pfq*(sJq.*fluxv)/2.0;

global useWADG MJ Mref
for e = 1:K
    if useWADG
        
        rp = Vq1D*reshape(rhsp(:,e),N+1,N+1)*Vq1D';
        ru = Vq1D*reshape(rhsu(:,e),N+1,N+1)*Vq1D';
        rv = Vq1D*reshape(rhsv(:,e),N+1,N+1)*Vq1D';
        JqK = reshape(Jq(:,e),Nq,Nq);
        
        rp = Pq1D*(rp./JqK)*Pq1D';
        ru = Pq1D*(ru./JqK)*Pq1D';
        rv = Pq1D*(rv./JqK)*Pq1D';
        rhsp(:,e) = rp(:);
        rhsu(:,e) = ru(:);
        rhsv(:,e) = rv(:);
    else
        rhsp(:,e) = MJ{e}\(Mref*rhsp(:,e));
        rhsu(:,e) = MJ{e}\(Mref*rhsu(:,e));
        rhsv(:,e) = MJ{e}\(Mref*rhsv(:,e));
        
    end
end

end
