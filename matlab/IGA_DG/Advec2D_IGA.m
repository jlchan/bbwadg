clear -globals
clear

Globals2D;
NB = 4;
Ksub = 32;
K1D = 1;
smoothKnots = 0;

global tau
tau = 1;

N = NB+Ksub-1;
dofs = (N+1)^2*K1D^2;

FinalTime = 2;

% Read in Mesh
[Nv, VX, VY, K, EToV] = QuadMesh2D(K1D);

a = 0;
ids = find(abs(abs(VX)-1) > 1e-8 & abs(abs(VY)-1) > 1e-8);
VX(ids) = VX(ids) + a;
VY(ids) = VY(ids) + a;

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
a = 0*.125;
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

%% build periodic maps

% make quadrature face maps

xf = x(Fmask(:),:);
yf = y(Fmask(:),:);

DX = max(VX)-min(VX);
DY = max(VY)-min(VY);

nbfaces = length(mapB)/(Nfp);
xB = reshape(xf(mapB),Nfp,nbfaces);
yB = reshape(yf(mapB),Nfp,nbfaces);
mapB = reshape(mapB,Nfp,nbfaces);
mapPtmp = mapP;
for i = 1:nbfaces
   for j = 1:nbfaces
       tol = 1e-8;
%        [x1 x2] = meshgrid(xB(:,i),xB(:,j));
%        [y1 y2] = meshgrid(yB(:,i),yB(:,j));
%        abs(x1-x2)
       dx = abs(sort(xB(:,i))-sort(xB(:,j)));
       dy = abs(sort(yB(:,i))-sort(yB(:,j)));
       
       if max(dx)<tol && max(abs(dy-DY))<tol
           mapPtmp(mapB(:,i)) = mapM(mapB(:,j));
%            keyboard
       elseif max(dy)<tol && max(abs(dx-DX))<tol
           mapPtmp(mapB(:,i)) = mapM(mapB(:,j));
       end      
   end
end
mapP = mapPtmp;

% for i = 1:Nfp*Nfaces*K
%    clf
%    plot(xf,yf,'.')
%    hold on
%    plot(xf(mapM(i)),yf(mapM(i)),'x','markersize',16)
%    plot(xf(mapP(i)),yf(mapP(i)),'o','markersize',16)
%    pause
% end
% return

%% check eigs

if 0 && Np*K < 3000
    U = zeros(Np*K,1);
    A = zeros(Np*K);
    for i = 1:Np*K
        U(i) = 1;
        u = reshape(U,Np,K);
        rhsu = advecRHS2Dq(u);
        A(:,i) = rhsu(:);
        U(i) = 0;
        if (mod(i,ceil(Np*K/10))==0)
            disp(sprintf('on i = %d out of %d\n',i,Np*K))
        end
    end
    lam = eig(A);
    plot(lam,'o');
    title(sprintf('largest real part = %g\n',max(real(lam))))
    max(abs(lam))
    keyboard
end

%% estimate timestep

U = randn((N+1)^2*K,1);
for i = 1:10
    Uprev = U;
    U = advecRHS2Dq(reshape(U,Np,K));
    
    lam = Uprev(:)'*U(:) / norm(Uprev(:));
    U = U/norm(U(:));
end
dt = .15/abs(lam);

dt = 10/(N*K1D*Ksub);

%% initial cond

x0 = 0; y0 = 0;
uex = @(x,y,t) exp(-5^2*((x-x0).^2 + (y-y0).^2));

u = zeros(Np,K);
for e = 1:K
    u(:,e) = reshape(Pq1D*reshape(uex(xq(:,e),yq(:,e),0),Nq,Nq)*Pq1D',(N+1)^2,1);
end

% vv = zeros(length(rp),K);
% for e = 1:K
%     vv(:,e) = reshape(Vp1D*reshape(u(:,e),N+1,N+1)*Vp1D',length(rp),1);
% end
% color_line3(xp,yp,vv,vv,'.');
% axis equal
% axis tight
% colorbar
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
        rhsu = advecRHS2Dq(u);

        % initiate and increment Runge-Kutta residuals
        resu = rk4a(INTRK)*resu + dt*rhsu;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        
    end;
    
    if (mod(tstep,1)==0 || tstep==Nsteps)
        clf
        vv = zeros(length(rp),K);
        for e = 1:K
            vv(:,e) = reshape(Vp1D*reshape(u(:,e),N+1,N+1)*Vp1D',length(rp),1);
        end
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f',tstep*dt))
        drawnow
    end
    
    if mod(tstep,25)==0
        disp(sprintf('on tstep %d out of %d\n',tstep,Nsteps))
    end
end


function rhs = advecRHS2Dq(u)

Globals2D;

global rxJ sxJ ryJ syJ nxq nyq Jq sJq
% global Vq Vrq Vsq Pq Prq Psq
global Vfq Pfq
global Vq1D Vrq1D Pq1D Prq1D

% Define field differences at faces
uf = u(Fmask(:),:);
du = zeros(Nfp*Nfaces,K); du(:) = uf(mapP)-uf(mapM);
du = Vfq*du;

bnu = nxq.*du;

global tau
flux =  tau.*(nxq.*nxq).*du - bnu;

Nq = size(Vrq1D,1);
ux = zeros(Np,K);
for e = 1:K
    urK = Vrq1D*reshape(u(:,e),N+1,N+1)*Vq1D';
    usK = Vq1D*reshape(u(:,e),N+1,N+1)*Vrq1D';        
    uxK = Pq1D*(reshape(rxJ(:,e),Nq,Nq).*urK + reshape(sxJ(:,e),Nq,Nq).*usK)*Pq1D';        
    ux(:,e) = uxK(:);
end


% compute right hand sides of the PDE's
rhs =  -ux + Pfq*(sJq.*flux)/2.0;

for e = 1:K
    rr = Vq1D*reshape(rhs(:,e),N+1,N+1)*Vq1D';
    JqK = reshape(Jq(:,e),Nq,Nq);
    
    rr = Pq1D*(rr./JqK)*Pq1D';
    rhs(:,e) = rr(:);
end

end
