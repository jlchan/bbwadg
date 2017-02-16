function [L2err] = WaveQuad_IGA(NB,Ksub,K1D,dt)

Globals2D;
if nargin==0
    NB = 4;
    Ksub = NB;
    K1D = 4;
    dt = .25/(NB*Ksub*K1D);
end

N = NB+Ksub-1

FinalTime = .7/pi;

% Read in Mesh
[Nv, VX, VY, K, EToV] = QuadMesh2D(K1D);

StartUpQuad2D;

% BuildPeriodicMaps2D(2,2);
% plot(x(vmapM),y(vmapM),'o')
% hold on
% plot(x(vmapP),y(vmapP),'x')

global rxJ sxJ ryJ syJ nxJ nyJ J 
global Vq Vfq Vrq Vsq Pq Prq Psq Pfq

a = .0;
x = x + a*cos(.5*3*pi*y).*(1-x).*(1+x);
y = y + a*sin(.5*3*pi*x).*(1-y).*(1+y);
[rx,sx,ry,sy,J] = GeometricFactorsQuad2D(x,y,Dr,Ds);
[nx, ny, sJ] = NormalsQuad2D();
Fscale = sJ./(J(Fmask,:));
% plot(x,y,'o');return

rp1D = linspace(-1,1,50);
[rp sp] = meshgrid(rp1D);
rp = rp(:); sp = sp(:);

Vp = Vandermonde2DQuad(N,rp,sp)/V;
xp = Vp*x;
yp = Vp*y;

%% initialize TP operators

[rq1D wq1D] = JacobiGQ(0,0,N+1);

[r1D w1D] = JacobiGL(0,0,N);
VX = linspace(-1,1,Ksub+1);
t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec
for i = 1:N+1
    r1D(i) = mean(t((i+1):(i+NB))); % greville
end

% switch to new nodal points (assumes boundary points still!)
[r s] = meshgrid(r1D); r = r(:); s = s(:);
Vnodal = Vandermonde1D(N,r1D)/Vandermonde1D(N,JacobiGL(0,0,N));
Vnodal = kron(Vnodal,Vnodal);
x = Vnodal*x;
y = Vnodal*y;

[rq sq] = meshgrid(rq1D); rq = rq(:); sq = sq(:);
[wrq wsq] = meshgrid(wq1D); wq = wrq(:).*wsq(:);
V1D = Vandermonde1D(N,r1D);
Vq1D = Vandermonde1D(N,rq1D)/V1D;
Vq = kron(Vq1D,Vq1D);
xq = Vq*x; yq = Vq*y;
wJq = diag(wq)*(Vq*J);

if 0
    V1D = Vandermonde1D(N,r1D);
    Dr1D = GradVandermonde1D(N,r1D)/V1D;
    M1D = inv(V1D*V1D');
    invM1D = V1D*V1D';
    
    M1D = diag(w1D);
    invM1D = diag(1./w1D);
    
    VDM = eye(Np);
    
else    
    
    [BVDM M1D Dr1D] = bsplineVDM(NB,Ksub,r1D); % VDM for interp, mass, M\S
    VDM = kron(BVDM,BVDM);
    Vp1D = bsplineVDM(NB,Ksub,rp1D);
    Vq1D = bsplineVDM(NB,Ksub,rq1D);
    
    % if collocating, change basis
    if 0
        M1D = inv(BVDM)'*M1D*inv(BVDM);
        Dr1D = BVDM*(Dr1D/BVDM);
        Vp1D = Vp1D/BVDM;
        Vq1D = Vq1D/BVDM;
        VDM = eye(Np);
        
        %         % optional: mass lump?
        %         M1D = diag(sum(M1D,2)); % diagonal norm
        %         M1D = diag(w1D);

    end
    invM1D = inv(M1D);
    Vp = kron(Vp1D,Vp1D);
    
end

Dr = kron(eye(N+1),Dr1D);
Ds = kron(Dr1D,eye(N+1));

M = kron(M1D,M1D);
invM = kron(invM1D,invM1D);
Mf = zeros(N+1,4*(N+1));
for f = 1:4
    Mf(Fmask(:,f),(1:N+1) + (f-1)*(N+1)) = M1D;
end
LIFT = invM * Mf;

global Vq Vfq Vrq Vsq Pq Prq Psq Pfq
Vq = kron(Vq1D,Vq1D);
Vfq = blkdiag(Vq1D,Vq1D,Vq1D,Vq1D); 
Vrq1D = Vq1D*Dr1D;
Vrq = kron(Vrq1D,Vq1D);
Vsq = kron(Vq1D,Vrq1D);

Pq = invM*Vq'*diag(wq);
Prq = invM*Vrq'*diag(wq);
Psq = invM*Vsq'*diag(wq);
Pfq = invM*Vfq'*diag(wfq);

rxJ = Vq*(rx.*J);
sxJ = Vq*(sx.*J);
ryJ = Vq*(ry.*J);
syJ = Vq*(sy.*J);
J = Vq*J;
nxJ = Vfq*(nx.*sJ);
nyJ = Vfq*(ny.*sJ);


%% check eigs

if 0 && nargin==0 && 3*Np*K < 1000
    U = zeros(Np*K,3);
    A = zeros(3*Np*K);
    for i = 1:3*Np*K
        U(i) = 1;
        p = reshape(U(:,1),Np,K);
        u = reshape(U(:,2),Np,K);
        v = reshape(U(:,3),Np,K);
        [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v);
        A(:,i) = [rhsp(:);rhsu(:);rhsv(:)];
        U(i) = 0;
        if (mod(i,ceil(3*Np*K/10))==0)
            disp(sprintf('on i = %d out of %d\n',i,3*Np*K))
        end
    end
    lam = eig(A);
    plot(lam,'o');
    title(sprintf('largest real part = %g\n',max(real(lam))))
    keyboard
end

%% initial cond

x0 = 0; y0 = .1;
pex = @(x,y,t) exp(-10^2*((x-x0).^2 + (y-y0).^2));
k = 7;
pex = @(x,y,t) cos(k*pi*x/2).*cos(k*pi*y/2).*cos(sqrt(2)*.5*k*pi*t);

p = VDM\pex(x,y,0);
u = zeros(Np, K); 
v = zeros(Np, K);

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
        
        [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v);
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        p = p+rk4b(INTRK)*resp;
        
    end;
    
    if 1 && nargin==0 && mod(tstep,10)==0
        clf
        vv = Vp*p;
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f',time))
        drawnow
    end
        
end


err = wJq.*(Vq*p - pex(xq,yq,FinalTime)).^2;
L2err = sqrt(sum(err(:)));


function [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

global rxJ sxJ ryJ syJ nxJ nyJ J 
global Vq Vfq Vrq Vsq Pq Prq Psq Pfq

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

uavg = zeros(Nfp*Nfaces,K); uavg(:) = u(vmapP)+u(vmapM);
vavg = zeros(Nfp*Nfaces,K); vavg(:) = v(vmapP)+v(vmapM);
ndotUavg = nx.*uavg + ny.*vavg;
ndotUavg(mapB) = 2*(u(vmapB).*nx(mapB) +  v(vmapB).*ny(mapB));

dp = Vfq*dp;
du = Vfq*du;
dv = Vfq*dv;

% Impose reflective boundary conditions (p+ = -p-)
ndotdU = nx.*du + ny.*dv;
ndotdU(mapB) = 0;
dp(mapB) = -2*p(vmapB);

tau = 1;
% fluxp =  tau*dp - ndotUavg;
fluxp =  tau*dp - ndotdU;
fluxu =  (tau*ndotdU - dp).*nx;
fluxv =  (tau*ndotdU - dp).*ny;

pr = Dr*p; ps = Ds*p;
dpdx = rxJ.*pr + sxJ.*ps;
dpdy = ryJ.*pr + syJ.*ps;
dudx = rxJ.*(Dr*u) + sxJ.*(Ds*u);
dvdy = ryJ.*(Dr*v) + syJ.*(Ds*v);
divU = dudx + dvdy;

% uq = Vq*u; vq = Vq*v;
% divU = -( Prq*(rxJ.*(uq) + ryJ.*(vq)) + Psq*(sxJ.*(uq) + syJ.*(vq)));

% compute right hand sides of the PDE's
rhsp =  -divU + LIFT*(sJ.*fluxp)/2.0;
rhsu =  -dpdx + LIFT*(sJ.*fluxu)/2.0;
rhsv =  -dpdy + LIFT*(sJ.*fluxv)/2.0;

rhsp = rhsp./J;
rhsu = rhsu./J;
rhsv = rhsv./J;

% rhsp = Pq*((Vq*rhsp)./J);
% rhsu = Pq*((Vq*rhsu)./J);
% rhsv = Pq*((Vq*rhsv)./J);
return;



