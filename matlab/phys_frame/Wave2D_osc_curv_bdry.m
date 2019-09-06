% clear all

Globals2D

N = 6; K1D = 10;
a = 0;
a = .05;
% a = .25;
CFL = .5;

useCurved = 1;
plotMesh = 0;

computeEigs = 0;
global tau
tau = 1;

FinalTime = 1;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
StartUp2D;

% BuildPeriodicMaps2D(2,0);
% PlotMesh2D; return

[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

%% make curvilinear mesh (still unstable?)

global Drq Dsq Vfqf Pfq Lq Prq Psq Vq Pq
global rxJ sxJ ryJ syJ nxq nyq sJq Jq Nfq  rxq sxq ryq syq


Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq)); % J's cancel out
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

Drq = Vq*Dr; Dsq = Vq*Ds;
[rq1D wq1D] = JacobiGQ(0,0,N);
V1D = Vandermonde1D(N,JacobiGL(0,0,N));
Vq1D = Vandermonde1D(N,rq1D)/V1D;
M1D = Vq1D'*diag(wq1D)*Vq1D;
Pq1D = M1D\(Vq1D'*diag(wq1D));
Pfq = kron(eye(Nfaces),Pq1D);
Vfqf = kron(eye(Nfaces),Vq1D);

rfq = [rq1D; -rq1D; -ones(size(rq1D))];
sfq = [-ones(size(rq1D)); rq1D; -rq1D];
% rfq = Vfqf*r(Fmask(:),:);
% sfq = Vfqf*s(Fmask(:),:);
wfq = [wq1D; wq1D; wq1D];

Vfq = Vandermonde2D(N,rfq,sfq)/V;
Lq = M\(Vfq'*diag(wfq)); % quadrature-based lift

Prq = M\(Drq'*diag(wq));
Psq = M\(Dsq'*diag(wq));

if useCurved
    
    kk = 25;
    w = .85;
    Hx = (tanh(kk*(x+w))-tanh(kk*(x-w)))/2;
    Hy = (tanh(kk*(y+w))-tanh(kk*(y-w)))/2;
    H = 1-Hx.*Hy;
%     color_line3(x,y,H,H,'.')
%     return

    % warp mesh
    kc = 3;
    Lx = 1; Ly = 1;
    x0 = 0; y0 = 0;
    x = x + Lx*a*cos(kc*1/2*pi*(x-x0)/Lx).*cos(kc*3/2*pi*(y-y0)/Ly).*H;
    y = y + Ly*a*sin(kc*3/2*pi*(x-x0)/Lx).*cos(kc*1/2*pi*(y-y0)/Ly).*H;    

    xq = Vq*x; yq = Vq*y;
    xp = Vp*x; yp = Vp*y;
    
    Nq = length(rq); Nfq = Nfp;
    [rxq,sxq,ryq,syq,Jq] = GeometricFactors2D(x,y,Vq*Dr,Vq*Ds);
    rxJ = rxq.*Jq;    sxJ = sxq.*Jq;
    ryJ = ryq.*Jq;    syJ = syq.*Jq;
    
    [rxf,sxf,ryf,syf,Jf] = GeometricFactors2D(x,y,Vfq*Dr,Vfq*Ds);
    rxJf = rxf.*Jf;    sxJf = sxf.*Jf;
    ryJf = ryf.*Jf;    syJf = syf.*Jf;
    
    nrJ = [-zeros(size(rq1D)); ones(size(rq1D)); -ones(size(rq1D)); ];
    nsJ = [-ones(size(rq1D)); ones(size(rq1D)); -zeros(size(rq1D)); ];
    nrJ = repmat(nrJ,1,K);
    nsJ = repmat(nsJ,1,K);
    
    nxJ = rxJf.*nrJ + sxJf.*nsJ;
    nyJ = ryJf.*nrJ + syJf.*nsJ;
    
    nxq = nxJ./Jf;
    nyq = nyJ./Jf;
    sJq = sqrt(nxq.^2 + nyq.^2);
    nxq = nxq./sJq; nyq = nyq./sJq;
    sJq = sJq.*Jf;
    
    if plotMesh        
        rp1D = linspace(-1,1,100)';
        Vp1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,JacobiGL(0,0,N));
        Vfp = kron(eye(Nfaces),Vp1D);
        xfp = Vfp*x(Fmask(:),:);
        yfp = Vfp*y(Fmask(:),:);
        figure
        plot(xfp,yfp,'k.')
        hold on
        xfq = Vfq*x;  yfq = Vfq*y;
        axis equal
        axis off
%         plot(xfq,yfq,'o'); hold on
%         quiver(xfq,yfq,nxq,nyq)
        L2err = nan;
        return
    end
    
end


% [nxq(:) nyq(:) sJq(:)],return
% testing

%% initial sol

k = 1; % frequency of solution
W = (2*k-1)/2*pi;
p = cos(W*x).*cos(W*y);

% x0 = 0; y0 = .25;
% p = exp(-200*((x-x0).^2 + (y-y0).^2));
u = zeros(Np, K); v = zeros(Np, K);

%% check eigenvalues of DG matrix
if computeEigs
    e = zeros(3*Np*K,1);
    A = zeros(3*Np*K);
    for i = 1:3*Np*K
        e(i) = 1;
        ids = 1:Np*K;
        p = reshape(e(ids),Np,K);
        u = reshape(e(ids + Np*K),Np,K);
        v = reshape(e(ids + 2*Np*K),Np,K);
        [rhsp, rhsu, rhsv] = acousticsRHS2D_skew(p,u,v,0);
        A(:,i) = [rhsp(:);rhsu(:);rhsv(:)];
        e(i) = 0;
        if mod(i,100)==0
            disp(sprintf('on column %d out of %d\n',i,3*Np*K))
        end
    end
    lam = eig(A);
    hold on;
    plot(lam,'o')
    title(sprintf('Largest real part = %e',max(real(lam))))
    
    return
    
    M = Vq'*diag(wq)*Vq; Mh = kron(spdiag(J(1,:)),M);
    
end


%%

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
dt = CFL*2/(CN*max(Fscale(:)));

% outer time step loop
tstep = 0;
figure
while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        timelocal = time + rk4c(INTRK)*dt;
        
        [rhsp, rhsu, rhsv] = acousticsRHS2D_skew(p,u,v,timelocal);
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        p = p+rk4b(INTRK)*resp;
        
    end;
    
    if mod(tstep,10)==0
        clf
        vv = Vp*p;
        pp = log(abs(vv));
        color_line3(xp,yp,vv,vv,'.');
        axis tight
        title(sprintf('Time = %g\n',time))
        drawnow
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
end

axis off

%%

%pex = @(x,y,time) cos(W*x).*cos(W*y).*cos(W*sqrt(2)*sqrt( )*time);
pex = @(x,y,time) cos(W*x).*cos(W*y).*cos(W*sqrt(2)*time);
pex = pex(xq,yq,FinalTime);
L2err2 = 0;
for e = 1:K
    diff = Vq*p(:,e) - pex(:,e);
    L2err2 = L2err2 + J(1,e)*sum(wq.*diff.^2);
end

L2err = sqrt(sum(L2err2))

% skew sym form
function [rhsp, rhsu, rhsv] = acousticsRHS2D_skew(p,u,v,timelocal)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

global Drq Dsq Vfqf Pfq Lq Prq Psq Vq Pq
global rxJ sxJ ryJ syJ nxq nyq sJq Jq Nfq rxq sxq ryq syq

Globals2D;

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);
uavg = zeros(Nfp*Nfaces,K); uavg(:) = .5*(u(vmapP) + u(vmapM));
vavg = zeros(Nfp*Nfaces,K); vavg(:) = .5*(v(vmapP) + v(vmapM));

% Impose reflective boundary conditions (p+ = -p-)
du(mapB) = 0; dv(mapB) = 0; dp(mapB) = -2*p(vmapB);


% can interp numerical fluxes *after* since mult by sJ is higher order
dp = Vfqf*dp;
du = Vfqf*du;
dv = Vfqf*dv;
uavg = Vfqf*uavg;
vavg = Vfqf*vavg;

% evaluate upwind fluxes
ndotdU = nxq.*du + nyq.*dv;
ndotavgU = nxq.*uavg + nyq.*vavg;

global tau
fluxp =  tau*dp - 2*ndotavgU;
fluxu =  (tau*ndotdU - dp).*nxq;
fluxv =  (tau*ndotdU - dp).*nyq;

% local derivatives of fields
pr = Drq*p; ps = Dsq*p;

dpdx = pr.*rxq + ps.*sxq;
dpdy = pr.*ryq + ps.*syq;

% compute right hand sides of the PDE's
%rhsp =  Pq*(-divU.*Jq) + Lq*(fluxp.*sJq/2.0);
uq = Vq*u;
vq = Vq*v;
Ur = rxq.*uq + ryq.*vq;
Us = sxq.*uq + syq.*vq;
rhsp =  Prq*(Ur.*Jq) + Psq*(Us.*Jq) + Lq*(fluxp.*sJq/2.0);
rhsu =  Pq*(-dpdx.*Jq) + Lq*(fluxu.*sJq/2.0);
rhsv =  Pq*(-dpdy.*Jq) + Lq*(fluxv.*sJq/2.0);

% apply inverse mass matrix
rhsp = Pq*((Vq*rhsp)./Jq);
rhsu = Pq*((Vq*rhsu)./Jq);
rhsv = Pq*((Vq*rhsv)./Jq);

end
