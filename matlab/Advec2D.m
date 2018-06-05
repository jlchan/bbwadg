function Advec2D


Globals2D

N = 4;
K1D = 8;
c_flag = 0;
FinalTime = 2;
cfun = @(x,y) ones(size(x));
% cfun = @(x,y) 1 + .5*sin(pi*x).*sin(pi*y); % smooth velocity
% cfun = @(x,y) (1 + .5*sin(2*pi*x).*sin(2*pi*y) + (y > 0)); % piecewise smooth velocity

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
StartUp2D;
BuildPeriodicMaps2D(2,2);

a = .075;
x = x + a*cos(pi/2*x).*cos(3*pi/2*y);
y = y + a*cos(3*pi/2*x).*cos(pi/2*y);
[rx sx ry sy J] = GeometricFactors2D(x,y,Dr,Ds);

if 1
    rp1D = linspace(-1,1,100)';
    Vp1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,JacobiGL(0,0,N));
    Vfp = kron(eye(Nfaces),Vp1D);
    xfp = Vfp*x(Fmask(:),:);
    yfp = Vfp*y(Fmask(:),:);
    plot(xfp,yfp,'k.')
    axis equal
    axis off
    
    return
end
% plot(x,y,'o');return

% plotting nodes
[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

global Vq Pq Pfq Pfqf Vfqf nxq nyq Fscaleq
Nq = 2*N+1;
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
%%
global bx by

bxex = @(x,y) ones(size(x));
byex = @(x,y) zeros(size(x));
bx = Pq*bxex(xq,yq);
by = Pq*byex(xq,yq);

%% params setup

x0 = 0; y0 = 0;

u = Pq*exp(-5^2*((xq-x0).^2 + (yq-y0).^2));


%%

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
CNh = max(CN*max(Fscale(:)));
dt = 2/CNh;

% outer time step loop
figure
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


function [rhsu] = advecRHS2D(u,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;
global bx by
global Vq Pq Pfq Pfqf Vfqf nxq nyq Fscaleq

tau = 1;
bxf = Vfqf*bx(Fmask(:),:);
byf = Vfqf*by(Fmask(:),:);
bn = bxf.*nxq + byf.*nyq;
du = Vfqf*reshape(u(vmapP)-u(vmapM),Nfp*Nfaces,K);
fluxu = Pfqf*(.5*du.*(tau*abs(bn)-bn));

fx = bx.*u; fy = by.*u;
dfx = rx.*(Dr*fx) + sx.*(Ds*fx);
dfy = ry.*(Dr*fy) + sy.*(Ds*fy);
divF = dfx + dfy;

%rhsu =  -divF + Pfq*(Fscaleq.*fluxu);
rhsu =  -divF + LIFT*(Fscale.*fluxu);


return;

