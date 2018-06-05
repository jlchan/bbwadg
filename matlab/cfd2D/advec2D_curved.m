function Advec2D


Globals2D

N = 4;
K1D = 16;
c_flag = 0;
FinalTime = .5;
CFL = 2;
a = 1/16;

[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/squareireg.neu');

% [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
StartUp2D;

BuildPeriodicMaps2D(2,2);

% plotting nodes
[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

global Vq Pq Pfq Pfqf Vfqf nxq nyq Fscaleq
Nq = 2*N;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq)); % J's cancel out
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

[rq1D wq1D] = JacobiGQ(0,0,N);
rfq = [rq1D; -rq1D; -ones(size(rq1D))];
sfq = [-ones(size(rq1D)); rq1D; -rq1D];
wfq = [wq1D; wq1D; wq1D];
Vq1D = Vandermonde1D(N,rq1D)/Vandermonde1D(N,JacobiGL(0,0,N));
% plot(rfq,sfq,'o')

Vfq = Vandermonde2D(N,rfq,sfq)/V;
Vfqf = kron(eye(3),Vq1D);
Mf = Vfq'*diag(wfq)*Vfq;
Pfq = M\(Vfq'*diag(wfq));

Pq1D = (Vq1D'*diag(wq1D)*Vq1D) \ (Vq1D'*diag(wq1D));
Pfqf = kron(eye(3),Pq1D);
 
nxq = Vfqf*nx;
nyq = Vfqf*ny;
Fscaleq = Vfqf*Fscale;

nrJ = [-zeros(size(rq1D)); ones(size(rq1D)); -ones(size(rq1D)); ];
nsJ = [-ones(size(rq1D)); ones(size(rq1D)); -zeros(size(rq1D)); ];



global rxJ sxJ ryJ syJ
rxJ = Vq*(rx.*J); sxJ = Vq*(sx.*J);
ryJ = Vq*(ry.*J); syJ = Vq*(sy.*J);
J = Vq*J;

sJ = Vfqf*sJ;

global Vrq Vsq
Vrq = Vq*Dr;
Vsq = Vq*Ds;

%% make curvilinear mesh (still unstable?)

x0 = 0; y0 = 0; 
Lx = 1; Ly = 1;
x = x + Lx*a*cos(1/2*pi*(x-x0)/Lx).*cos(3/2*pi*(y-y0)/Ly);
y = y + Ly*a*sin(3/2*pi*(x-x0)/Lx).*cos(1/2*pi*(y-y0)/Ly);

xq = Vq*x; yq = Vq*y;
xp = Vp*x; yp = Vp*y;

if 1
    rp1D = linspace(-1,1,100)';
    Vp1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,JacobiGL(0,0,N));
    Vfp = kron(eye(Nfaces),Vp1D);
    xfp = Vfp*x(Fmask(:),:);
    yfp = Vfp*y(Fmask(:),:);
    figure(1)
    plot(xfp,yfp,'k.')
%     return
end

Nq = length(rq);
Nfq = length(rfq)/Nfaces;


rxJ = zeros(Nq,K); sxJ = zeros(Nq,K);
ryJ = zeros(Nq,K); syJ = zeros(Nq,K);
J = zeros(Nq,K);
rxJf = zeros(Nfq*Nfaces,K); sxJf = zeros(Nfq*Nfaces,K);
ryJf = zeros(Nfq*Nfaces,K); syJf = zeros(Nfq*Nfaces,K);
Jf = zeros(Nfq*Nfaces,K);
for e = 1:K
    [rxk,sxk,ryk,syk,Jk] = GeometricFactors2D(x(:,e),y(:,e),Vq*Dr,Vq*Ds);
    rxJ(:,e) = rxk.*Jk;    sxJ(:,e) = sxk.*Jk;
    ryJ(:,e) = ryk.*Jk;    syJ(:,e) = syk.*Jk;
    J(:,e) = Jk;
    
    [rxk,sxk,ryk,syk,Jk] = GeometricFactors2D(x(:,e),y(:,e),Vfq*Dr,Vfq*Ds);
    rxJf(:,e) = rxk.*Jk;    sxJf(:,e) = sxk.*Jk;
    ryJf(:,e) = ryk.*Jk;    syJf(:,e) = syk.*Jk;
    Jf(:,e) = Jk;
end

nxJ = rxJf.*nrJ + sxJf.*nsJ;
nyJ = ryJf.*nrJ + syJf.*nsJ;

nx = nxJ./Jf;
ny = nyJ./Jf;
sJ = sqrt(nx.^2 + ny.^2);
nx = nx./sJ; ny = ny./sJ;
sJ = sJ.*Jf;
Fscale = sJ./Jf;

xf = Vfq*x;
yf = Vfq*y;

%%
global bx by

bxex = @(x,y) ones(size(x));
byex = @(x,y) zeros(size(x));
bx = Pq*bxex(xq,yq);
by = Pq*byex(xq,yq);

%% params setup

x0 = -.5; y0 = 0;

u = Pq*exp(-5^2*((xq-x0).^2 + (yq-y0).^2));

%%
if 1
    u = zeros(Np,K);
    for i = 1:Np*K
        u(i) = 1;
        rhsu = advecRHS2D(u,0);
        A(:,i) = rhsu(:);
        u(i) = 0;
    end
    M = kron(diag(J(1,:)),inv(V*V'));
%     A = M*A;
    A(abs(A)<1e-8) = 0;
    spy(A)
    figure;
    spy(M)
    return
end

%%

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
CNh = max(CN*max(Fscale(:)));
dt = CFL/CNh;

wJq = diag(wq)*J;

% outer time step loop
figure(2)
% colormap(gray)
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;
for i = 1:Nsteps    
    for INTRK = 1:5
        
        timelocal = time + rk4c(INTRK)*dt;
        
        rhsu  = advecRHS2D(u,timelocal);
        
        % initiate and increment Runge-Kutta residuals
        resu = rk4a(INTRK)*resu + dt*rhsu;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        
    end;
    
    l2norm(i) = sum(sum(((Vq*u).^2).*wJq));
    
    if mod(i,10)==0 || i==Nsteps
        clf
        pp = u;
        vv = Vp*pp;        
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f',dt*i))
        drawnow
    end
        
end

figure(3)
plot((1:Nsteps)*dt,l2norm,'--')
hold on


function [rhsu] = advecRHS2D(u,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;
global bx by
global Vq Pq Pfq Pfqf Vfqf nxq nyq Fscaleq
global rxJ sxJ ryJ syJ sJ
global Vrq Vsq 

tau = 1;
bxf = Vfqf*bx(Fmask(:),:);
byf = Vfqf*by(Fmask(:),:);
bn = bxf.*nxq + byf.*nyq;
du = Vfqf*reshape(u(vmapP)-u(vmapM),Nfp*Nfaces,K);
fluxu = .5*du.*(tau*abs(bn)-bn);

fx = bx.*u; fy = by.*u;
dfx = rxJ.*(Vrq*fx) + sxJ.*(Vsq*fx);
dfy = ryJ.*(Vrq*fy) + syJ.*(Vsq*fy);
divF = dfx + dfy;

rhsu =  -Pq*(divF./J) + Pfq*(Fscale.*fluxu);


return;

