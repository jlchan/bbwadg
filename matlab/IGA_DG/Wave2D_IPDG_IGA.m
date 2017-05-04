function [L2err] = WaveQuad_IGA(NBin,Ksubin,K1Din)

Globals2D;

global NB Ksub

if nargin==0
    NB = 4;
    Ksub = 4;
    K1D = 4;
else
    NB = NBin;
    Ksub = Ksubin;
    K1D = K1Din;
end

smoothKnots = 25;
useQuadrature = 1;
    
N = NB+Ksub-1;
dofs = (N+1)^2*K1D^2

FinalTime = .3;

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
global Vq Vfq Vrq Vsq Pq Prq Psq Pfq

% non-affine mappings
a = .125;
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
Vp = Vandermonde2DQuad(N,rp,sp)/V;
xp = Vp*x;
yp = Vp*y;

[rq1D wq1D] = JacobiGQ(0,0,N+1);
if (Ksub>1)
    %[rq1D wq1D] = spline_quadrature(NB);
    [rgq wgq] = JacobiGQ(0,0,NB+1);
    h = (2/Ksub);
    rq1D = (1+repmat(rgq,1,Ksub))/2*h + repmat(linspace(-1,1-h,Ksub),length(rgq),1);    
    wq1D = repmat(wgq,Ksub,1)*h/2;    
    rq1D = rq1D(:); wq1D = wq1D(:);
end

[r1D] = JacobiGL(0,0,N);
if Ksub > 1
    [~, ~, ~, ~, ~, ~, ~, ~, VX] = bsplineVDM(NB,Ksub,r1D,smoothKnots);
%     VX = linspace(-1,1,Ksub+1);
    t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec
    for i = 1:N+1
        r1D(i) = mean(t((i+1):(i+NB))); % greville
    end
end


% switch to new nodal points (assumes boundary points still!)
[r s] = meshgrid(r1D); r = r(:); s = s(:);
Vnodal1D = Vandermonde1D(N,r1D)/Vandermonde1D(N,JacobiGL(0,0,N));
Vnodal = kron(Vnodal1D,Vnodal1D);
Vfnodal = blkdiag(Vnodal1D,Vnodal1D,Vnodal1D,Vnodal1D);
x = Vnodal*x;
y = Vnodal*y;
J = Vnodal*J;
rx = Vnodal*rx;
ry = Vnodal*ry;
sx = Vnodal*sx;
sy = Vnodal*sy;
nx = Vfnodal * nx;
ny = Vfnodal * ny;
sJ = Vfnodal * sJ;


[rq sq] = meshgrid(rq1D); rq = rq(:); sq = sq(:);
[wrq wsq] = meshgrid(wq1D); wq = wrq(:).*wsq(:);
V1D = Vandermonde1D(N,r1D);
Vq1D = Vandermonde1D(N,rq1D)/V1D;
Vq = kron(Vq1D,Vq1D);
xq = Vq*x; 
yq = Vq*y;
wJq = diag(wq)*(Vq*J);
if useQuadrature
    Vfq = blkdiag(Vq1D,Vq1D,Vq1D,Vq1D);
    rxJ = Vq*(rx.*J); sxJ = Vq*(sx.*J);
    ryJ = Vq*(ry.*J); syJ = Vq*(sy.*J);
    Jq  = Vq*J;
    nxq = Vfq*nx; nyq = Vfq*ny;
    sJq = Vfq*sJ;    

    % vv = J(:,1);
    % color_line3(x(:,1),y(:,1),vv,vv,'.');
    % hold on
    % vv = Jq(:,1);
    % color_line3(xq(:,1),yq(:,1),vv,vv,'o');
    % return
end

if Ksub==1
    V1D = Vandermonde1D(N,r1D);
    D1D = GradVandermonde1D(N,r1D)/V1D;
    Vp1D = Vandermonde1D(N,rp1D)/V1D;
    Vq1D = Vandermonde1D(N,rq1D)/V1D;
    M1D = Vq1D'*diag(wq1D)*Vq1D;
    invM1D = V1D*V1D';
else
    [BVDM M1D D1D] = bsplineVDM(NB,Ksub,r1D,smoothKnots); % VDM for interp, mass, M\S    
    Vp1D = bsplineVDM(NB,Ksub,rp1D,smoothKnots);
    Vq1D = bsplineVDM(NB,Ksub,rq1D,smoothKnots);
    invM1D = inv(M1D);
end
Vp = kron(Vp1D,Vp1D);

Dr = kron(eye(N+1),D1D);
Ds = kron(D1D,eye(N+1));

invM = kron(invM1D,invM1D);
Mf = zeros(N+1,4*(N+1));
for f = 1:4
    Mf(Fmask(:,f),(1:N+1) + (f-1)*(N+1)) = M1D;
end
LIFT = invM * Mf;
LIFT(abs(LIFT)<1e-8) = 0; 
% spy(LIFT);return

global Vq Vfq Vrq Vsq Pq Prq Psq Pfq

Vq = kron(Vq1D,Vq1D);

if useQuadrature
    Vfq = blkdiag(Vq1D,Vq1D,Vq1D,Vq1D);
    Vrq1D = Vq1D*D1D;
    Vrq = kron(Vq1D,Vrq1D);
    Vsq = kron(Vrq1D,Vq1D);
    
    Pq = invM*Vq'*diag(wq);
    Prq = invM*Vrq'*diag(wq);
    Psq = invM*Vsq'*diag(wq);    
    
    % project onto face space then apply usual lift    
    Pfq = LIFT*kron(eye(Nfaces),invM1D)*Vfq'*diag(repmat(wq1D,Nfaces,1));
    
else
    rxJ = (rx.*J);
    sxJ = (sx.*J);
    ryJ = (ry.*J);
    syJ = (sy.*J);
    nxq = (nx.*sJ);
    nyq = (ny.*sJ);
end

global Grr Grs Gss Gst 
Grr = (rxJ.^2 + ryJ.^2)./(Jq.^2);
Grs = (rxJ.*sxJ + ryJ.*syJ)./(Jq.^2);
Gss = (sxJ.^2 + syJ.^2)./(Jq.^2);

% keyboard

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
for e = 1:K
    p(:,e) = (Vq'*diag(wq.*Jq(:,e))*Vq)\(Vq'* (wq.*Jq(:,e).*pex(xq(:,e),yq(:,e),0)));
end

u = zeros(Np, K); 

err = wJq.*(Vq*p - pex(xq,yq,0)).^2;
init_cond_err = sqrt(sum(err(:)))
% return
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
        if useQuadrature
            [rhsu] = acousticsRHS2Dq(p);
            rhsp = u;
        else
            [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v);
        end        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        p = p+rk4b(INTRK)*resp;        
        
    end;
    
    if 1 && nargin==0 && (mod(tstep,5)==0 || tstep==Nsteps)
        clf
        vv = Vp*p;
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


err = wJq.*(Vq*p - pex(xq,yq,FinalTime)).^2;
L2err = sqrt(sum(err(:)));


function [rhs] = acousticsRHS2Dq(p)

Globals2D;

global rxJ sxJ ryJ syJ nxq nyq Jq sJq
global Vq Vfq Vrq Vsq Pq Prq Psq Pfq

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapM) - p(vmapP);
dp(mapB) = 2*p(vmapB);
dpq = Vfq*dp;

% compute q
prq = Vrq*p; 
psq = Vsq*p; 

%ux = Pq*((Vq*Pq*(rxJ.*prq + sxJ.*psq))./Jq);
%uy = Pq*((Vq*Pq*(ryJ.*prq + syJ.*psq))./Jq);
ux = Pq*((rxJ.*prq + sxJ.*psq)./Jq);
uy = Pq*((ryJ.*prq + syJ.*psq)./Jq);

qx = ux - .5*Pfq*(dpq.*nxq.*sJq);
qy = uy - .5*Pfq*(dpq.*nyq.*sJq);

dqx = Vfq*reshape(qx(vmapM)-.5*(ux(vmapM)+ux(vmapP)),Nfp*Nfaces,K);
dqy = Vfq*reshape(qy(vmapM)-.5*(uy(vmapM)+uy(vmapP)),Nfp*Nfaces,K);

qxr = Vrq*qx; qxs = Vsq*qx;
qyr = Vrq*qy; qys = Vsq*qy;

global NB Ksub
CT = max((NB+1)*Ksub,(NB+1)^2);
hmin = min(min(Jq)./max(sJq));
tau = CT / hmin;
dqn = nxq.*dqx + nyq.*dqy;
fluxq = dqn + tau*dpq;

rhs = Pq*(rxJ.*qxr + sxJ.*qxs + ryJ.*qyr + syJ.*qys) - Pfq*(sJq.*fluxq);

rhs = Pq*((Vq*rhs)./Jq);

