function ElasticAcoustic2D_unified

% clear all, clear
clear -global *

Globals2D

K1D = 64;
N = 5;
c_flag = 0;
FinalTime = .5;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

StartUp2D;

[rp sp] = EquiNodes2D(15); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
Pq = V*V'*Vq'*diag(wq); % J's cancel out
Mref = inv(V*V');
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

% partition mesh: y > 0 = acoustic
global Ka Ke
% Ke = 1:K; Ka = []; 
% Ka = 1:K; Ke = [];
Ka = find(mean(y) > 0); Ke = find(mean(y) < 0);

%% find dividing boundary

global mapAM vmapAM mapAP vmapAP
global mapEM vmapEM mapEP vmapEP
mapAM = [];    mapAP = [];
vmapAM = [];    vmapAP = [];
mapEM = [];    mapEP = [];
vmapEM = [];    vmapEP = [];

if ~isempty(Ke) && ~isempty(Ka)
    
    mapAM = zeros(Nfp,K1D);    mapAP = zeros(Nfp,K1D);
    vmapAM = zeros(Nfp,K1D);    vmapAP = zeros(Nfp,K1D);
    mapEM = zeros(Nfp,K1D);    mapEP = zeros(Nfp,K1D);
    vmapEM = zeros(Nfp,K1D);    vmapEP = zeros(Nfp,K1D);
    
    yM = reshape(y(vmapM),Nfp*Nfaces,K);
    ska = 1; ske = 1;
    for e = 1:K
        yf = reshape(yM(:,e),Nfp,Nfaces);
        for f = 1:Nfaces
            if norm(yf(:,f))<1e-8
                ids = (1:Nfp) + (f-1)*Nfp + (e-1)*Nfp*Nfaces;
                if ismember(e,Ka) % if in acoustic region
                    mapAM(:,ska) = mapM(ids);
                    mapAP(:,ska) = mapP(ids);
                    vmapAM(:,ska) = vmapM(ids);
                    vmapAP(:,ska) = vmapP(ids);
                    ska = ska + 1;
                elseif ismember(e,Ke) % if in elastic region
                    mapEM(:,ske) = mapM(ids);
                    mapEP(:,ske) = mapP(ids);
                    vmapEM(:,ske) = vmapM(ids);
                    vmapEP(:,ske) = vmapP(ids);
                    ske = ske + 1;
                else
                    keyboard
                end                
            end
        end
    end
    
%     keyboard
end

%%
global Nfld mu lambda Vq Pq tau c2
Nfld = 5; %(u1,u2,sxx,syy,sxy)

mu = ones(size(xq));
lambda = ones(size(xq));
c2 = ones(size(xq));

k = 3;
mu = 1 + .5*cos(k*pi*xq).*cos(k*pi*yq);
c2 = 1 + .5*cos(k*pi*xq).*cos(k*pi*yq);

mu = V\(Pq*mu); 
c2 = V\(Pq*c2);
mu = repmat(mu(1,:),length(wq),1)/sqrt(2);
c2 = repmat(c2(1,:),length(wq),1)/sqrt(2);

% vv = Vp*Pq*c2; color_line3(xp,yp,vv,vv,'.');return
keyboard

tau = 1;

%% params setup

x0 = 0; 
y0 = .1;
pp = exp(-100^2*((x-x0).^2 + (y-y0).^2));

f0 = 10;
t0 = 1/f0;

global fsrc
fsrc = @(t) (t < t0).*(1-2*(pi*f0*(t-t0))^2)*exp(-(pi*f0*(t-t0)^2)).* (Pq * exp(-100^2*((xq-x0).^2 + (yq-y0).^2)));
% y0 = -.25;
% p = p + exp(-10^2*((x-x0).^2 + (y-y0).^2));

u = zeros(Np, K);

U{1} = u;
U{2} = u;
U{3} = u;
U{4} = u;
U{5} = u;

%%
if Nfld*Np*K < 2500
        
    u = zeros(Nfld*Np*K,1);
    rhs = zeros(Nfld*Np*K,1);
    A = zeros(Nfld*Np*K);
    ids = 1:Np*K;
    for i = 1:Nfld*Np*K
        u(i) = 1;
        for fld = 1:Nfld
            U{fld} = reshape(u(ids + (fld-1)*Np*K),Np,K);
        end
        rE = ElasRHS2D(U,0);
        rA = AcousRHS2D(U,0);        
        
        u(i) = 0;
        for fld = 1:Nfld  
            rU = zeros(Np,K);
            rU(:,Ke) = rE{fld}(:,Ke);
            if fld <= 3
                rU(:,Ka) = rA{fld}(:,Ka);
            end
            rhs(ids + (fld-1)*Np*K) = rU;
        end
        
        A(:,i) = rhs(:);
        if (mod(i,100)==0)
            disp(sprintf('on col %d out of %d\n',i,Np*K*Nfld))
        end
    end
    lam = eig(A);
    
    plot(lam,'.','markersize',24)
    hold on
    title(sprintf('Largest real part = %g\n',max(real(lam))))
    axis equal
    %         drawnow
    max(abs(lam))
    
    keyboard
end
%%

time = 0;

% Runge-Kutta residual storage
for fld = 1:Nfld
    res{fld} = zeros(Np,K);
end

% compute time step size
CN = (N+1)*(N+2)/3; % guessing...
dt = 2/(max(2*mu(:)+lambda(:))*CN*max(Fscale(:)));
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

% outer time step loop
tstep = 0;

M = inv(V*V');
wqJ = diag(wq)*(Vq*J);

% figure
% colormap(gray)
% colormap(hot)
for tstep = 1:Nsteps
        
    time = tstep*dt;
    
    for INTRK = 1:5
        
        timeloc = time + rk4c(INTRK)*dt;
        rhsE = ElasRHS2D(U,timeloc);
        rhsA = AcousRHS2D(U,timeloc);        
        
        % initiate and increment Runge-Kutta residuals
        for fld = 1:Nfld
            rhs = zeros(Np,K);                    
            if fld <= 3
                rhs(:,Ka) = rhsA{fld}(:,Ka);                
            end            
            rhs(:,Ke) = rhsE{fld}(:,Ke);
            res{fld} = rk4a(INTRK)*res{fld} + dt*rhs;
            U{fld} = U{fld} + rk4b(INTRK)*res{fld};
        end
        
    end;
    
    if 1 && (mod(tstep,10)==0 || tstep==Nsteps)
        clf
        
        pe = (U{3} + U{4})/2; % trace(S)
% pe = U{3}; % trace(S)
        pa = U{3}; % trace(S)
        p(:,Ke) = pe(:,Ke);
        p(:,Ka) = pa(:,Ka);
        vv = Vp*p;
        color_line3(xp,yp,vv,vv,'.');
        axis tight
        title(sprintf('time = %f',time));
        colorbar;
        
        drawnow        
        
    end
    
    if mod(tstep,100)==0
        disp(sprintf('On timestep %d out of %d\n',tstep,round(FinalTime/dt)))
    end
end

% keyboard
set(gca,'fontsize',14)
title('')
axis tight


function [rhs] = ElasRHS2D(U,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

global Nfld mu lambda Vq Pq tau useWADG
global C11 C12 C13 C22 C23 C33

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

% stress fluxes
nUx = dU{1}.*nx;
nUy = dU{2}.*ny;
nUxy = dU{2}.*nx + dU{1}.*ny;

opt=1;
if opt==1 % traction BCs    
    nSx(mapB) = -2*(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -2*(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));    
elseif opt==2 % basic ABCs
    nSx(mapB) = -(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
    dU{1}(mapB) = -U{1}(vmapB);
    dU{2}(mapB) = -U{2}(vmapB);
elseif opt==3 % zero velocity
    dU{1}(mapB) = -2*U{1}(vmapB);
    dU{2}(mapB) = -2*U{2}(vmapB);    
end

% coupling
global mapEM vmapEM mapEP vmapEP 
nxf = nx(mapEM); nyf = ny(mapEM); 
u = U{1}(vmapEP); 
v = U{2}(vmapEP);
p = U{3}(vmapEP);

Snx = U{3}(vmapEM).*nxf + U{5}(vmapEM).*nyf;
Sny = U{5}(vmapEM).*nxf + U{4}(vmapEM).*nyf;
% dU1 = (u-U{1}(vmapEM));
% dU2 = (v-U{2}(vmapEM));
% dUn = dU1.*nxf + dU2.*nyf;
% dU1 = dUn;
% dU2 = dUn;
Un = u.*nxf + v.*nyf;
dU1 = (Un.*nxf-U{1}(vmapEM));
dU2 = (Un.*nyf-U{2}(vmapEM));

nSx(mapEM) = (p.*nxf - Snx);
nSy(mapEM) = (p.*nyf - Sny);
nUx(mapEM) = dU1.*nxf;
nUy(mapEM) = dU2.*nyf;
nUxy(mapEM) = dU1.*nyf + dU2.*nxf;

% nSx(mapEM) = -2*Snx;
% nSy(mapEM) = -2*Sny;
% nUx(mapEM) = 0;
% nUy(mapEM) = 0;
% nUxy(mapEM) = 0;

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
    flux{fld}(:) = fc{fld}(:) + tau.*fp{fld}(:);
end

% compute right hand sides of the PDE's
rr{1} =  divSx   +  LIFT*(Fscale.*flux{1})/2.0;
rr{2} =  divSy   +  LIFT*(Fscale.*flux{2})/2.0;
rr{3} =  du1dx   +  LIFT*(Fscale.*flux{3})/2.0;
rr{4} =  du2dy   +  LIFT*(Fscale.*flux{4})/2.0;
rr{5} =  du12dxy +  LIFT*(Fscale.*flux{5})/2.0;

if 0
    rhs{1} = rr{1};
    rhs{2} = rr{2};
    rhs{3} = (2*mu+lambda).*rr{3} + lambda.*rr{4};
    rhs{4} = lambda.*rr{3} + (2*mu+lambda).*rr{4};
    rhs{5} = (mu) .* rr{5};
else
    global Pq Vq
    rhs{1} = rr{1};
    rhs{2} = rr{2};
    rhs{3} = Pq*((2*mu+lambda).*(Vq*rr{3}) + lambda.*(Vq*rr{4}));
    rhs{4} = Pq*(lambda.*(Vq*rr{3}) + (2*mu+lambda).*(Vq*rr{4}));
    rhs{5} = Pq*(mu .* (Vq*rr{5}));
end

% global Ka
% for fld = 1:5
%     rhs{fld}(:,Ka) = 0;
% end

return;


function [rhs] = AcousRHS2D(U,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

u = U{1};
v = U{2};
p = U{3}; % should be same as U{4}

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

% evaluate upwind fluxes
ndotdU = nx.*du + ny.*dv;

% Impose reflective boundary conditions (p+ = -p-)
ndotdU(mapB) = 0;
dp(mapB) = -2*p(vmapB);

% elastic-acoustic coupling
global mapAM vmapAM mapAP vmapAP 
nxf = nx(mapAM); nyf = ny(mapAM); 
v1 = U{1}(vmapAP); v2 = U{2}(vmapAP);
sxx = U{3}(vmapAP); syy = U{4}(vmapAP); sxy = U{5}(vmapAP);

Snx = sxx.*nxf + sxy.*nyf;
Sny = sxy.*nxf + syy.*nyf;
nSn = Snx.*nxf + Sny.*nyf;
dU1 = (v1-u(vmapAM));
dU2 = (v2-v(vmapAM));
ndotdU(mapAM) = dU1.*nxf + dU2.*nyf;
dp(mapAM) = nSn-p(vmapAM);

% ndotdU(mapAM) = 0;
% dp(mapAM) = -2*p(vmapAM);

global tau;
fluxp =  tau*dp + ndotdU;
fluxu =  (tau*ndotdU + dp).*nx;
fluxv =  (tau*ndotdU + dp).*ny;

pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
rhs{1} = dpdx + LIFT*(Fscale.*fluxu)/2.0;
rhs{2} = dpdy + LIFT*(Fscale.*fluxv)/2.0;
rhs{3} = divU + LIFT*(Fscale.*fluxp)/2.0;

global c2 Pq Vq

global fsrc 
rhs{3} = rhs{3} + fsrc(time);

rhs{3} = Pq*(c2.*(Vq*rhs{3}));

