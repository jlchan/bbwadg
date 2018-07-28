% function L2err = Elasticity2D_Lamb(N,K1D)

% clear all, clear
clear -global *
clear
Globals2D

K1D = 2;
N = 3;
FinalTime = 1;
tau0 = 0;
CFL = 1;

a = .0;
useCurved = 1;
computeEigs = 1;

% filename = 'Grid/Other/block2.neu';
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
[Nv, VX, VY, K, EToV] = unif_tri_mesh(2*K1D,K1D);
VY = (VY+1)/2-.5;

StartUp2D;

BuildPeriodicMaps2D(2,1);
% BuildPeriodicMaps2D(2,0);

[rp sp] = EquiNodes2D(15); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq)); % J's cancel out
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;


%% make curvilinear mesh (still unstable?)

global Drq Dsq Vfqf Pfq Lq Prq Psq
global rxJ sxJ ryJ syJ nxq nyq sJq Jq Nfq

Drq = Vq*Dr; Dsq = Vq*Ds;
[rq1D wq1D] = JacobiGQ(0,0,N);
V1D = Vandermonde1D(N,JacobiGL(0,0,N));
Vq1D = Vandermonde1D(N,rq1D)/V1D;
M1D = Vq1D'*diag(wq1D)*Vq1D;
Pq1D = M1D\(Vq1D'*diag(wq1D));
Pfq = kron(eye(Nfaces),Pq1D);
Vfqf = kron(eye(Nfaces),Vq1D);

% rfq = [rq1D; -rq1D; -ones(size(rq1D))];
% sfq = [-ones(size(rq1D)); rq1D; -rq1D];
rfq = Vfqf*r(Fmask(:),:);
sfq = Vfqf*s(Fmask(:),:);
wfq = [wq1D; wq1D; wq1D];
%wfq = [wq1D; wq1D; wq1D*sqrt(2)];
% wfq = wfq/sum(wfq);

Vfq = Vandermonde2D(N,rfq,sfq)/V;
Lq = M\(Vfq'*diag(wfq)); % quadrature-based lift

Prq = M\(Drq'*diag(wq));
Psq = M\(Dsq'*diag(wq));

if useCurved
    
    % warp mesh
    Lx = 1; Ly = 1/2;
    x0 = 0; y0 = 0;
    x = x + Lx*a*cos(1/2*pi*(x-x0)/Lx).*cos(3/2*pi*(y-y0)/Ly);
    y = y + Ly*a*sin(3/2*pi*(x-x0)/Lx).*cos(1/2*pi*(y-y0)/Ly);
    
    xq = Vq*x; yq = Vq*y;
    xp = Vp*x; yp = Vp*y;
    
    Nq = length(rq); Nfq = Nfp;
    [rxq,sxq,ryq,syq,Jq] = GeometricFactors2D(x,y,Vq*Dr,Vq*Ds);
    rxJ = rxq.*Jq;    sxJ = sxq.*Jq;
    ryJ = ryq.*Jq;    syJ = syq.*Jq;
    [rxf,sxf,ryf,syf,Jf] = GeometricFactors2D(x,y,Vfq*Dr,Vfq*Ds);
    rxJf = rxf.*Jf;    sxJf = sxf.*Jf;
    ryJf = ryf.*Jf;    syJf = syf.*Jf;
    
    if 1
        nrJ = [-zeros(size(rq1D)); ones(size(rq1D)); -ones(size(rq1D)); ];
        nsJ = [-ones(size(rq1D)); ones(size(rq1D)); -zeros(size(rq1D)); ];
%         plot(rfq,sfq,'o'); hold on; quiver(rfq,sfq,nrJ,nsJ); return
        nxJ = rxJf.*nrJ + sxJf.*nsJ;
        nyJ = ryJf.*nrJ + syJf.*nsJ;
        
        nxq = nxJ./Jf;
        nyq = nyJ./Jf;
        sJq = sqrt(nxq.^2 + nyq.^2);
        nxq = nxq./sJq; nyq = nyq./sJq;
        sJq = sJq.*Jf;
    end
    
    if 0
        rp1D = linspace(-1,1,100)';
        Vp1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,JacobiGL(0,0,N));
        Vfp = kron(eye(Nfaces),Vp1D);
        xfp = Vfp*x(Fmask(:),:);
        yfp = Vfp*y(Fmask(:),:);
        plot(xfp,yfp,'k.')
        hold on
        xfq = Vfq*x;  yfq = Vfq*y;
        plot(xfq,yfq,'o'); hold on
        quiver(xfq,yfq,nx,ny)
        L2err = nan;
        return
    end
    
end

% testing
vmapP = reshape(vmapP,Nfp*Nfaces,K);
vmapM = reshape(vmapM,Nfp*Nfaces,K);

p = randn(Np,K);
u = randn(Np,K);
v = randn(Np,K);
pf = Vfqf*p(Fmask(:),:);
uf = Vfqf*u(Fmask(:),:);
vf = Vfqf*v(Fmask(:),:);
% Vfqf = eye(size(Vfqf));
pjump = Vfqf*(p(vmapP)-p(vmapM));

% pf = p(Fmask(:),:);
% pjump2 = Vfqf*reshape(pf(mapP)-pf(mapM),Nfp*Nfaces,K);
% return

ujump = Vfqf*(u(vmapP)-u(vmapM));
uavg = .5*Vfqf*(u(vmapP)+u(vmapM));

vjump = Vfqf*(v(vmapP)-v(vmapM));
vavg = .5*Vfqf*(v(vmapP)+v(vmapM));

unf = uf.*nxq + vf.*nyq;
unavg = uavg.*nxq + vavg.*nyq;

wsJ = diag(wfq)*sJq;
Wf = diag(wsJ(:));

% (nxq(:).*uf(:) + nyq(:).*vf(:))'*Wf*(pf(mapP)-pf(mapM))
% unf(:)'*Wf*pjump(:)
% unavg(:)'*Wf*pjump(:)
% 2*pf(:)'*Wf*unavg(:)

cf_error = abs(unf(:)'*Wf*pjump(:) + 2*pf(:)'*Wf*unavg(:)) % central flux error

return
    


%% problem params

global Nfld mu lambda rho tau useWADG
Nfld = 5; %(u1,u2,sxx,syy,sxy)

rho = 1;
lambda = 2;
mu = 1;

for fld = 1:Nfld
    tau{fld} = tau0;
end


%% check eigs


if computeEigs
    
    u = zeros(Nfld*Np*K,1);
    rhs = zeros(Nfld*Np*K,1);
    A = zeros(Nfld*Np*K);
    ids = 1:Np*K;
    for i = 1:Nfld*Np*K
        u(i) = 1;
        for fld = 1:Nfld
            U{fld} = reshape(u(ids + (fld-1)*Np*K),Np,K);
        end
        if useCurved
            rU = ElasRHS2D(U,0);            
        else
            rU = ElasRHS2D_affine(U,0);
        end
        u(i) = 0;
        for fld = 1:Nfld
            rhs(ids + (fld-1)*Np*K) = rU{fld};
        end
        A(:,i) = rhs(:);
        if (mod(i,100)==0)
            disp(sprintf('on col %d out of %d\n',i,Np*K*Nfld))
        end
    end
    lam = eig(A);
    plot(lam,'.','markersize',16)
    hold on
    title(sprintf('Largest real part = %g\n',max(real(lam))))
%     axis equal
    %         drawnow
    max(abs(lam))
    
    return
end

%% params setup

% cp = sqrt((2*mu + lambda)/rho);
% cs = sqrt(mu/rho);

% assume constant
global k B1 p q B2 w

k = 2*pi;
w = 13.137063197233;
B1 = 126.1992721468;
B2 = 53.88807700007;
p = sqrt(w^2/(2*mu+lambda) - k^2);
q = sqrt(w^2/mu - k^2);

u1 = @(x,y,t) (-k*B1*cos(p*y) - q*B2*cos(q*y)).*sin(k*x-w*t);
u2 = @(x,y,t) (-p*B1*sin(p*y) + k*B2*sin(q*y)).*cos(k*x-w*t);

v1 = @(x,y,t) (B1*k*cos(p*y) + B2*q*cos(q*y)).*(w*cos(k*x - t*w));
v2 = @(x,y,t) (B2*k*sin(q*y) - B1*p*sin(p*y)).*(w*sin(k*x - t*w));

u1x = @(x,y,t) -k.*cos(k.*x - t.*w).*(B1.*k.*cos(p.*y) + B2.*q.*cos(q.*y));
u2y = @(x,y,t) -cos(k.*x - t.*w).*(B1.*p^2.*cos(p.*y) - B2.*k.*q.*cos(q.*y));
u12xy = @(x,y,t) sin(k.*x - t.*w).*(B2.*q^2.*sin(q.*y) - B2.*k^2.*sin(q.*y) + 2.*B1.*k.*p.*sin(p.*y));

sxx = @(x,y,t) ((2*mu+lambda) .* u1x(x,y,t)) + (lambda.*u2y(x,y,t));
syy = @(x,y,t) lambda.*u1x(x,y,t) + ((2*mu+lambda) .* u2y(x,y,t));
sxy = @(x,y,t) mu .* u12xy(x,y,t);

U{1} = v1(x,y,0);
U{2} = v2(x,y,0);
U{3} = ((2*mu + lambda)*u1x(x,y,0) + lambda*u2y(x,y,0));
U{4} = (lambda*u1x(x,y,0) + (2*mu + lambda)*u2y(x,y,0));
U{5} = mu*u12xy(x,y,0);

useWADG = 1;
if useWADG
    mu = mu*ones(size(xq));
    lambda = lambda*ones(size(xq));
end

%%

time = 0;

% Runge-Kutta residual storage
for fld = 1:Nfld
    res{fld} = zeros(Np,K);
end

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
h = 1/K1D;
dt = CFL*h/(max(mu(:)+lambda(:))*CN);

% outer time step loop
tstep = 0;


while (time<FinalTime)
    
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        timeloc = time + rk4c(INTRK)*dt;
        if useCurved
            rhs = ElasRHS2D(U,timeloc);
            
%             if tstep==10
%                 val = zeros(5,1);
%                 for fld = 1:5
%                     val(fld) = sum(sum(U{fld}.*(M*rhs{fld})));
%                 end                
%                 keyboard
%             end

        else
            rhs = ElasRHS2D_affine(U,timeloc);
        end
        
        % initiate and increment Runge-Kutta residuals
        for fld = 1:Nfld
            res{fld} = rk4a(INTRK)*res{fld} + dt*rhs{fld};
            U{fld} = U{fld} + rk4b(INTRK)*res{fld};
        end
        
    end;
    
    if mod(tstep,10)==0
        clf
        
        %         p = U{3}+U{4}; % trace(S)
        vp1 = Vp*U{1};
        vp1 = vp1/2.5e4;
        color_line3(xp,yp,vp1,vp1,'.');
        axis([-1 1 -.5 .5 -1.5 1.5])
        
        title(sprintf('time = %f',time))
        colorbar
        drawnow
        
        
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
    if mod(tstep,100)==0
        disp(sprintf('On timestep %d out of %d\n',tstep,round(FinalTime/dt)))
    end
end


Nq = 3*N;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

wqJ = diag(wq)*Jq;
err = (Vq*U{1} - v1(xq,yq,FinalTime)).^2 + (Vq*U{2} - v2(xq,yq,FinalTime)).^2;
uq = v1(xq,yq,FinalTime).^2 + v2(xq,yq,FinalTime).^2;
errU = sqrt(sum(wqJ(:).*err(:)))/sqrt(sum(wqJ(:).*uq(:)));

serr = (Vq*U{3} - sxx(xq,yq,FinalTime)).^2 + (Vq*U{4} - syy(xq,yq,FinalTime)).^2 + (Vq*U{5} - sxy(xq,yq,FinalTime)).^2;
sq = sxx(xq,yq,FinalTime).^2 + syy(xq,yq,FinalTime).^2 + sxy(xq,yq,FinalTime).^2;
errS = sqrt(sum(wqJ(:).*serr(:)))/sqrt(sum(wqJ(:).*sq(:)));

L2err = sqrt(sum(wqJ(:).*(err(:)+serr(:))))/sqrt(sum(wqJ(:).*(uq(:)+sq(:))))
%  title(sprintf('L2 err = %g\n',L2err))

% vv = Vp*p;
% err = (abs(u1(xp,yp,FinalTime)-vv));
% figure
% color_line3(xp,yp,err,err,'.')




function [rhs] = ElasRHS2D_affine(U,time)

Globals2D;

global Nfld mu lambda rho Vq Pq tau useWADG
global mapBd vmapBd mapBt vmapBt

% Define field differences at faces
for fld = 1:Nfld
    u = U{fld};
    
    % compute jumps
    dU{fld} = zeros(Nfp*Nfaces,K);
    dU{fld}(:) = u(vmapP)-u(vmapM);
    
    ur = Dr*u;
    us = Ds*u;
    Ux{fld} = rx.*ur + sx.*us;
    Uy{fld} = ry.*ur + sy.*us;
end

divSx = Ux{3} + Uy{5}; % d(Sxx)dx + d(Sxy)dy
divSy = Ux{5} + Uy{4}; % d(Sxy)dx + d(Syy)dy
du1dx = Ux{1}; % du1dx
du2dy = Uy{2}; % du2dy
du12dxy = Ux{2} + Uy{1}; % du2dx + du1dy

% velocity fluxes
nSx = nx.*dU{3} + ny.*dU{5};
nSy = nx.*dU{5} + ny.*dU{4};

% traction BCs
nSx(mapB) = -2*(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
nSy(mapB) = -2*(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));

% stress fluxes
nUx = dU{1}.*nx;
nUy = dU{2}.*ny;
nUxy = dU{2}.*nx + dU{1}.*ny;

% evaluate central fluxes
fc{1} = nSx;
fc{2} = nSy;
fc{3} = nUx;
fc{4} = nUy;
fc{5} = nUxy;

% penalization terms - reapply An
fp{1} = nx.*fc{3} + ny.*fc{5};
fp{2} = nx.*fc{5} + ny.*fc{4};
fp{3} = fc{1}.*nx;
fp{4} = fc{2}.*ny;
fp{5} = fc{2}.*nx + fc{1}.*ny;

flux = cell(5,1);
for fld = 1:Nfld
    flux{fld} = zeros(Nfp*Nfaces,K);
    flux{fld}(:) = fc{fld}(:) + tau{fld}.*fp{fld}(:);
end

% compute right hand sides of the PDE's
rr{1} =  divSx   +  LIFT*(Fscale.*flux{1})/2.0;
rr{2} =  divSy   +  LIFT*(Fscale.*flux{2})/2.0;
rr{3} =  du1dx   +  LIFT*(Fscale.*flux{3})/2.0;
rr{4} =  du2dy   +  LIFT*(Fscale.*flux{4})/2.0;
rr{5} =  du12dxy +  LIFT*(Fscale.*flux{5})/2.0;

% C = [2*mu+lambda       lambda       0
%      lambda       2*mu+lambda       0
%           0       0                mu]
if useWADG
    for fld = 3:Nfld
        rr{fld} = Vq*rr{fld};
    end
    rhs{1} = (1/rho)*rr{1};
    rhs{2} = (1/rho)*rr{2};
    rhs{3} = Pq*((2*mu+lambda).*rr{3} + lambda.*rr{4});
    rhs{4} = Pq*(lambda.*rr{3} + (2*mu+lambda).*rr{4});
    rhs{5} = Pq*((mu) .* rr{5});
else
    rhs{1} = (1/rho)*rr{1};
    rhs{2} = (1/rho)*rr{2};
    rhs{3} = (2*mu+lambda).*rr{3} + lambda.*rr{4};
    rhs{4} = lambda.*rr{3} + (2*mu+lambda).*rr{4};
    rhs{5} = (mu) .* rr{5};
end

end


function [rhs] = ElasRHS2D(U,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

global Nfld mu lambda rho tau useWADG
global mapBd vmapBd mapBt vmapBt
global Drq Dsq Vfqf Pfq Lq Prq Psq Vq Pq
global rxJ sxJ ryJ syJ nxq nyq sJq Jq Nfq

% Define field differences at faces
for fld = 1:Nfld
    u = U{fld};
    
    dU{fld} = zeros(Nfp*Nfaces,K);
    dU{fld}(:) = Vfqf*reshape(u(vmapP)-u(vmapM),Nfp*Nfaces,K);
    
    if fld <= 2
        aU{fld} = zeros(Nfp*Nfaces,K);
        aU{fld}(:) = Vfqf*reshape(u(vmapP) + u(vmapM),Nfp*Nfaces,K);
        uq = Vq*u;
        Ux{fld} = -(Prq*(rxJ.*uq) + Psq*(sxJ.*uq));
        Uy{fld} = -(Prq*(ryJ.*uq) + Psq*(syJ.*uq));
    else
        aU{fld} = zeros(Nfp*Nfaces,K);
        aU{fld}(:) = Vfqf*reshape(u(vmapP)-u(vmapM),Nfp*Nfaces,K);
        
        ur = Drq*u;  us = Dsq*u;
        Ux{fld} = Pq*(rxJ.*ur + sxJ.*us);
        Uy{fld} = Pq*(ryJ.*ur + syJ.*us);
    end
end

divSx = Ux{3} + Uy{5}; % d(Sxx)dx + d(Sxy)dy
divSy = Ux{5} + Uy{4}; % d(Sxy)dx + d(Syy)dy
du1dx = Ux{1}; % du1dx
du2dy = Uy{2}; % du2dy
du12dxy = Ux{2} + Uy{1}; % du2dx + du1dy

% velocity fluxes
nSx = nxq.*dU{3} + nyq.*dU{5};
nSy = nxq.*dU{5} + nyq.*dU{4};

% traction BCs
U3M = Vfqf*reshape(U{3}(vmapM),Nfp*Nfaces,K);
U4M = Vfqf*reshape(U{4}(vmapM),Nfp*Nfaces,K);
U5M = Vfqf*reshape(U{5}(vmapM),Nfp*Nfaces,K);
U3B = U3M(mapB); U4B = U4M(mapB); U5B = U5M(mapB);
nSx(mapB) = -2*(nxq(mapB).*U3B + nyq(mapB).*U5B);
nSy(mapB) = -2*(nxq(mapB).*U5B + nyq(mapB).*U4B);

% evaluate central fluxes
fc{1} = nSx;
fc{2} = nSy;
fc{3} = aU{1}.*nxq;
fc{4} = aU{2}.*nyq;
fc{5} = aU{2}.*nxq + aU{1}.*nyq;

% penalization terms
nUxy = dU{2}.*nxq + dU{1}.*nyq;
fp{1} = nxq.^2.*dU{1}  + nyq.*nUxy;
fp{2} = nxq.*nUxy + nyq.^2.*dU{2};
fp{3} = nSx.*nxq;
fp{4} = nSy.*nyq;
fp{5} = nSy.*nxq + nSx.*nyq;

flux = cell(5,1);
for fld = 1:Nfld
    flux{fld} = zeros(Nfp*Nfaces,K);
    flux{fld}(:) = fc{fld}(:) + tau{fld}.*fp{fld}(:);
end

% compute right hand sides of the PDE's
rr{1} =  divSx   +  Lq*(sJq.*flux{1})/2.0;
rr{2} =  divSy   +  Lq*(sJq.*flux{2})/2.0;
rr{3} =  du1dx   +  Lq*(sJq.*flux{3})/2.0;
rr{4} =  du2dy   +  Lq*(sJq.*flux{4})/2.0;
rr{5} =  du12dxy +  Lq*(sJq.*flux{5})/2.0;

% C = [2*mu+lambda       lambda       0
%      lambda       2*mu+lambda       0
%           0       0                mu]
for fld = 1:Nfld
    rhs{fld} = rr{fld};
end

% for fld = 1:Nfld
%     rr{fld} = Vq*rr{fld};
% end
% rhs{1} = Pq*(rr{1}./(rho.*Jq));
% rhs{2} = Pq*(rr{2}./(rho.*Jq));
% rhs{3} = Pq*(((2*mu+lambda).*rr{3} + lambda.*rr{4})./Jq);
% rhs{4} = Pq*((lambda.*rr{3} + (2*mu+lambda).*rr{4})./Jq);
% rhs{5} = Pq*((mu .* rr{5})./Jq);

end
