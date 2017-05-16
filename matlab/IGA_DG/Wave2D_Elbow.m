function [L2err] = WaveQuad_IGA(NB,Ksub,K1D,dt)

Globals2D;
if nargin==0
    NB = 5;
    Ksub = 16;
    K1D = 1;
end

smoothKnots = 0;
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
global Vq1D Vrq1D Pq1D Prq1D Vfq

% initialize TP operators
rp1D = linspace(-1,1,150);
[rp sp] = meshgrid(rp1D);
rp = rp(:); sp = sp(:);

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

[rq sq] = meshgrid(rq1D); rq = rq(:); sq = sq(:);
[wrq wsq] = meshgrid(wq1D); wq = wrq(:).*wsq(:);
V1D = Vandermonde1D(N,r1D);
Vq1D = Vandermonde1D(N,rq1D)/V1D;
% Vq = kron(Vq1D,Vq1D);

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
% Vp = kron(Vp1D,Vp1D);

% Dr = kron(eye(N+1),D1D);
% Ds = kron(D1D,eye(N+1));

invM = kron(invM1D,invM1D);
Mf = zeros(N+1,4*(N+1));
for f = 1:4
    Mf(Fmask(:,f),(1:N+1) + (f-1)*(N+1)) = M1D;
end
LIFT = invM * Mf;
LIFT(abs(LIFT)<1e-8) = 0; 
LIFT = sparse(LIFT);
% spy(LIFT);return
% keyboard

% global Prq Psq 
if useQuadrature
    Vrq1D = Vq1D*D1D;
%     Vfq = blkdiag(Vq1D,Vq1D,Vq1D,Vq1D);            
    Pq1D = invM1D*Vq1D'*diag(wq1D);
    Prq1D = invM1D*Vrq1D'*diag(wq1D);
       
%     Vq = kron(Vq1D,Vq1D);
%     Pq = invM*Vq'*diag(wq);
%     Vrq = kron(Vq1D,Vrq1D);
%     Vsq = kron(Vrq1D,Vq1D);
%     Prq = invM*Vrq'*diag(wq);
%     Psq = invM*Vsq'*diag(wq);        
%     % project onto face space then apply usual lift    
%     Pfq = LIFT*kron(eye(Nfaces),invM1D)*Vfq'*diag(repmat(wq1D,Nfaces,1));   
   
end

disp(sprintf('making geom ...\n'))
[xp yp] = Elbow2D(rp1D);
[xq yq rxq sxq ryq syq Jq] = Elbow2D(rq1D);

e = ones(size(rq1D));
rfq = [rq1D e rq1D -e];
sfq = [-e rq1D e rq1D];

xfq = []; yfq = [];
rxf = []; sxf = []; ryf = []; syf = []; 
Jf = [];
for f = 1:4
    if f==1
        [xff yff rxff sxff ryff syff Jff] = Elbow2D(rq1D,-1);
    elseif f==2
        [xff yff rxff sxff ryff syff Jff] = Elbow2D(1,rq1D);
    elseif f==3
        [xff yff rxff sxff ryff syff Jff] = Elbow2D(rq1D,1);
    elseif f==4
        [xff yff rxff sxff ryff syff Jff] = Elbow2D(-1,rq1D);
    end
    xfq = [xfq xff];
    yfq = [yfq yff];
    rxf = [rxf rxff];
    sxf = [sxf sxff];
    ryf = [ryf ryff];
    syf = [syf syff];
    Jf = [Jf Jff];
end

disp('done making geom')

nxq = zeros(size(rxf));
nyq = zeros(size(rxf));

nxq(:,1) = -sxf(:,1); nyq(:,1) = -syf(:,1);
nxq(:,2) =  rxf(:,2); nyq(:,2) = ryf(:,2);
nxq(:,3) =  sxf(:,3); nyq(:,3) = syf(:,3);
nxq(:,4) = -rxf(:,4); nyq(:,4) = -ryf(:,4);

sJq = sqrt(nxq.*nxq+nyq.*nyq);
nxq = nxq./sJq; nyq = nyq./sJq;
sJq = sJq.*Jf;

% % plot(rq,sq,'o')
% hold on
% % plot(rfq,sfq,'o')
% % keyboard
% plot(xfq,yfq,'o'); hold on; quiver(xfq,yfq,nxq,nyq) ;return
% color_line3(xq,yq,Jq,Jq,'.'); return

nxq = nxq(:); nyq = nyq(:); sJq = sJq(:);
rxJ = rxq.*Jq; sxJ = sxq.*Jq;
ryJ = ryq.*Jq; syJ = syq.*Jq;
wJq = spdiag(wq)*Jq;


%% check eigs

3*Np*K
if 0 && nargin==0 && 3*Np*K < 3000
    U = zeros(Np*K,3);
    A = zeros(3*Np*K);
    for i = 1:3*Np*K
        U(i) = 1;
        p = reshape(U(:,1),Np,K);
        u = reshape(U(:,2),Np,K);
        v = reshape(U(:,3),Np,K);
        if useQuadrature
            [rhsp, rhsu, rhsv] = acousticsRHS2Dq(p,u,v);
        else
            [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v);
        end
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

%% estimate timestep
if 1
    U = randn((N+1)^2*K,3);
    for i = 1:10
        Uprev = U;
        if useQuadrature
            [rhsp, rhsu, rhsv] = acousticsRHS2Dq(reshape(U(:,1),Np,K),reshape(U(:,2),Np,K),reshape(U(:,3),Np,K));
        else
            [rhsp, rhsu, rhsv] = acousticsRHS2D(reshape(U(:,1),Np,K),reshape(U(:,2),Np,K),reshape(U(:,3),Np,K));
        end
        U(:,1) = rhsp(:);
        U(:,2) = rhsu(:);
        U(:,3) = rhsv(:);
        
        lam = Uprev(:)'*U(:) / norm(Uprev(:));
        U = U/norm(U(:));
    end
    dt = .75/abs(lam)
else
    CT = max((NB+1)*Ksub,(NB+1)^2);
    dt = .25/(CT*K1D);
end

%% initial cond

x0 = -.75*cos(pi/4); y0 = .75*sin(pi/4);
% x0 = 0; y0 = 0;
pex = @(x,y,t) exp(-10^2*((x-x0).^2 + (y-y0).^2));
% k = 3;
% pex = @(x,y,t) cos(k*pi*x/2).*cos(k*pi*y/2).*cos(sqrt(2)*.5*k*pi*t);

Nq = size(Pq1D,2);
p = Pq1D*reshape(pex(xq,yq,0),Nq,Nq)*Pq1D'; p  = p(:);
% for e = 1:K
%     p(:,e) = (Vq'*diag(wq.*Jq(:,e))*Vq)\(Vq'* (wq.*Jq(:,e).*pex(xq(:,e),yq(:,e),0)));
% end

% p = VDM\pex(x,y,0);
u = zeros(Np, K); 
v = zeros(Np, K);

pq = Vq1D*reshape(p,N+1,N+1)*Vq1D'; pq = pq(:);
err = wJq.*(pq - pex(xq,yq,0)).^2;
init_cond_err = sqrt(sum(err(:)))
% return

% vv = Vp1D*reshape(p,N+1,N+1)*Vp1D';
% color_line3(xp,yp,vv,vv,'.')
% return
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
        if useQuadrature
            [rhsp, rhsu, rhsv] = acousticsRHS2Dq(p,u,v);
        else
            [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v);
        end        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        
        % update fields
        p = p+rk4b(INTRK)*resp;
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        
        
    end;
    
    if 1 && nargin==0 && (mod(tstep,5)==0 || tstep==Nsteps)
        clf
        vv = Vp1D*reshape(p,N+1,N+1)*Vp1D';
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

return

% err = wJq.*(Vq*p - pex(xq,yq,FinalTime)).^2;
% L2err = sqrt(sum(err(:)));


function [rhsp, rhsu, rhsv] = acousticsRHS2Dq(p,u,v)

Globals2D;

global rxJ sxJ ryJ syJ nxq nyq Jq sJq

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

uavg = zeros(Nfp*Nfaces,K); uavg(:) = u(vmapP)+u(vmapM);
vavg = zeros(Nfp*Nfaces,K); vavg(:) = v(vmapP)+v(vmapM);

% uavg(mapB) = u(vmapB);
% vavg(mapB) = v(vmapB);
% Impose reflective boundary conditions (p+ = -p-)
% dp(mapB) = -2*p(vmapB);
du(mapB) = -2*u(vmapB);
dv(mapB) = -2*v(vmapB);
uavg(mapB) = 0;
vavg(mapB) = 0;

global Vq1D Vrq1D Pq1D Prq1D 
Nq = size(Vq1D,1);
dp = Vq1D*reshape(dp,N+1,4); dp = dp(:);
du = Vq1D*reshape(du,N+1,4); du = du(:);
dv = Vq1D*reshape(dv,N+1,4); dv = dv(:);
uavg = Vq1D*reshape(uavg,N+1,4); uavg = uavg(:);
vavg = Vq1D*reshape(vavg,N+1,4); vavg = vavg(:);


ndotUavg = nxq.*uavg + nyq.*vavg;
ndotdU = nxq.*du + nyq.*dv;

tau = 1;
fluxp =  tau*dp - ndotUavg;
fluxu =  (tau*ndotdU - dp).*nxq;
fluxv =  (tau*ndotdU - dp).*nyq;

% uq = Vq*u; vq = Vq*v;
uq = Vq1D*reshape(u,N+1,N+1)*Vq1D'; uq = uq(:);
vq = Vq1D*reshape(v,N+1,N+1)*Vq1D'; vq = vq(:);

Ur = rxJ.*uq + ryJ.*vq;
Us = sxJ.*uq + syJ.*vq;
PUr = (Prq1D*reshape(Ur,Nq,Nq))*Pq1D';
PUs = (Pq1D*reshape(Us,Nq,Nq))*Prq1D';
divU = -(PUr + PUs); divU = divU(:);

%rhsp =  (-divU) + Pfq*(sJq.*fluxp)/2.0;
fp = Pq1D*reshape(sJq.*fluxp,Nq,4); fp = fp(:);
rhsp =  (-divU) + .5*LIFT*(fp);

pr = Vrq1D*reshape(p,N+1,N+1)*Vq1D'; pr = pr(:);
ps = Vq1D*reshape(p,N+1,N+1)*Vrq1D'; ps = ps(:);

dpdx = rxJ.*pr + sxJ.*ps;
dpdy = ryJ.*pr + syJ.*ps;
px = Pq1D*reshape(dpdx,Nq,Nq)*Pq1D'; px = px(:);
py = Pq1D*reshape(dpdy,Nq,Nq)*Pq1D'; py = py(:);

fu = Pq1D*reshape(sJq.*fluxu,Nq,4); fu = fu(:);
fv = Pq1D*reshape(sJq.*fluxv,Nq,4); fv = fv(:);
rhsu = -px + .5*LIFT*(fu);
rhsv = -py + .5*LIFT*(fv);
%rhsu = -px + Pfq*(sJq.*fluxu)/2.0;
%rhsv = -py + Pfq*(sJq.*fluxv)/2.0;

rhsp = Pq1D*((Vq1D*reshape(rhsp,N+1,N+1)*Vq1D') ./ reshape(Jq,Nq,Nq))*Pq1D'; rhsp = rhsp(:);
rhsu = Pq1D*((Vq1D*reshape(rhsu,N+1,N+1)*Vq1D') ./ reshape(Jq,Nq,Nq))*Pq1D'; rhsu = rhsu(:);
rhsv = Pq1D*((Vq1D*reshape(rhsv,N+1,N+1)*Vq1D') ./ reshape(Jq,Nq,Nq))*Pq1D'; rhsv = rhsv(:);

