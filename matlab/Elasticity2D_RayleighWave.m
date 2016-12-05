function Elasticity2D_RayleighWave

% clear all, clear
clear -global *

Globals2D

mu = .1;
lambda = 1;
FinalTime = .25*sqrt(1/mu);

K1D = 4;
N = 3;

% filename = 'Grid/Other/block2.neu';
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
aa = 2;
% [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D,aa*K1D);
% VX = (VX+1)/2; VY = -(VY+1)/2*aa;
[Nv, VX, VY, K, EToV] = unif_tri_mesh(aa*K1D,K1D);
VX = aa*(VX+1)/2; VY = (VY+1)/2;

StartUp2D;

% PlotMesh2D;return

BuildPeriodicMaps2D(0,1);

% for i = 1:length(mapB)
%     clf
%     plot(x,y,'.')
%     hold on
%     id = vmapM(mapB(i));
%    plot(x(id), y(id),'o','markersize',32)
%    id = vmapP(mapB(i));
%    plot(x(id), y(id),'x','markersize',32)
%    pause
% end
% return

[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
Pq = V*V'*Vq'*diag(wq); % J's cancel out
Mref = inv(V*V');
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

%% BC stuff

global mapBd vmapBd mapBt vmapBt

% find x = 0 faces
fbids = reshape(find(vmapP==vmapM),Nfp,nnz(vmapP==vmapM)/Nfp);
% yfb = y(vmapM(fbids)); bfaces = sum(abs(yfb),1)<1e-8;
xfb = x(vmapM(fbids)); bfaces = sum(abs(xfb),1)<1e-8;

fbids = fbids(:,bfaces);
fbids = fbids(:);

mapBt = fbids;
vmapBt = vmapP(fbids);

% remove Dirichlet BCs for traction
mapBd = setdiff(mapB,mapBt);
vmapBd = vmapP(mapBd);

if 0
    plot(x(vmapBt),y(vmapBt),'bo')
    hold on
    plot(x(vmapBd),y(vmapBd),'rd')
    plot(x,y,'k.')
    return
end
%%
global Nfld mu lambda Vq Pq tau useWADG
Nfld = 5; 

rho = 1;

% C = [2*mu+lambda       lambda       0
%      lambda       2*mu+lambda       0
%           0       0                mu];

% plot3(x,y,mu,'.');return
tau0 = 1;
for fld = 1:5
    tau{fld} = tau0;
    if fld > 2
        tau{fld} = tau0./(mu+lambda);
    end
end

%% params setup

% assume constant
global v1 v2

% format long
d = @(xi) sqrt(1-xi.^2).*sqrt(1-xi.^2*mu / (2*mu + lambda)) - (.5*xi.^2-1).^2;
xi = fzero(d,.9,optimset('TolX',1e-12));

cr = xi*sqrt(mu);
w = 2*pi;

a = @(x) exp(-w*x*sqrt(1-xi^2));
b = @(x) (.5*xi^2 - 1)*exp(-w*x*sqrt(1-xi^2*mu / (2*mu+lambda)));

u1 = @(x,y,t) (a(x) + b(x)).*cos(w*(y+cr.*t));
u2 = @(x,y,t) (a(x).*sqrt(1-xi^2) + b(x)./ sqrt(1-xi^2*mu / (2*mu+lambda))).*sin(w*(y+cr*t));

u1x = @(x,y,t) -cos(w.*(y + cr.*t)).*(w.*exp(-w.*x.*(1 - xi^2)^(1/2)).*(1 - xi^2)^(1/2) + w.*exp(-w.*x.*(1 - (mu.*xi^2)/(lambda + 2.*mu))^(1/2)).*(xi^2/2 - 1).*(1 - (mu.*xi^2)/(lambda + 2.*mu))^(1/2));
u2y = @(x,y,t) w.*cos(w.*(y + cr.*t)).*(exp(-w.*x.*(1 - xi^2)^(1/2)).*(1 - xi^2)^(1/2) + (exp(-w.*x.*(1 - (mu.*xi^2)/(lambda + 2.*mu))^(1/2)).*(xi^2/2 - 1))/(1 - (mu.*xi^2)/(lambda + 2.*mu))^(1/2));
u12xy = @(x,y,t)-sin(w.*(y + cr.*t)).*(w.*exp(-w.*x.*(1 - (mu.*xi^2)/(lambda + 2.*mu))^(1/2)).*(xi^2/2 - 1) - w.*exp(-w.*x.*(1 - xi^2)^(1/2)).*(xi^2 - 1)) - w.*sin(w.*(y + cr.*t)).*(exp(-w.*x.*(1 - xi^2)^(1/2)) + exp(-w.*x.*(1 - (mu.*xi^2)/(lambda + 2.*mu))^(1/2)).*(xi^2/2 - 1));

v1 = @(x,y,t) -(a(x)+b(x)).*(cr*w*sin(w*(y + cr*t)));
v2 = @(x,y,t) (a(x).*sqrt(1-xi^2) + b(x)./ sqrt(1-xi^2*mu / (2*mu+lambda))).*(cr*w*cos(w*(y + cr*t)));


sxx = @(x,y,t) ((2*mu+lambda) .* u1x(x,y,t)) + (lambda.*u2y(x,y,t));
syy = @(x,y,t) lambda.*u1x(x,y,t) + ((2*mu+lambda) .* u2y(x,y,t));
sxy = @(x,y,t) mu .* u12xy(x,y,t);

U{1} = v1(x,y,0);
U{2} = v2(x,y,0);
U{3} = (2*mu+lambda) * u1x(x,y,0) + lambda * u2y(x,y,0);
U{4} = lambda * u1x(x,y,0) + (2*mu+lambda) * u2y(x,y,0);
U{5} = mu * (u12xy(x,y,0));

% for fld = 3:5
%     U{fld} = zeros(Np,K);
% end

if 0
    for t = 0:.1:4
        clf
        %         subplot(2,1,1)
        v1 = sxx(xp,yp,t)+ syy(xp,yp,t);
        v1 = v1/10;
        color_line3(xp,yp,v1,v1,'.')
        axis equal
        view(3)
        %         view(90,0)
        
        %         subplot(2,1,2)
        %
        %         color_line3(xp+v1,yp+v2,v2,v2,'.')
        %         axis equal
        %         title(sprintf('time = %f\n',t))
        drawnow
        
    end
    return
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
dt = 1/(max(2*mu(:)+lambda(:))*CN*max(Fscale(:)));

% outer time step loop
tstep = 0;


figure
% colormap(gray)
% colormap(hot)
while (time<FinalTime)
    
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        timeloc = time + rk4c(INTRK)*dt;
        rhs = ElasRHS2D(U,timeloc);
        
        % initiate and increment Runge-Kutta residuals
        for fld = 1:Nfld
            res{fld} = rk4a(INTRK)*res{fld} + dt*rhs{fld};
            U{fld} = U{fld} + rk4b(INTRK)*res{fld};
        end
        
    end;
    
    if 1 && mod(tstep,10)==0
        clf
%                 subplot(1,2,1)
%         vv = Vp*U{1}/10;
        vv = Vp*(U{3}+U{4})/10;
        color_line3(xp,yp,vv,vv,'.');
        axis tight
        title(sprintf('time = %f',time))
        colorbar
        axis([0 aa 0 1 -1 1])

        
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

wqJ = diag(wq)*(Vq*J);
err = (Vq*U{1} - v1(xq,yq,FinalTime)).^2 + (Vq*U{2} - v2(xq,yq,FinalTime)).^2;
uq = v1(xq,yq,FinalTime).^2 + v2(xq,yq,FinalTime).^2;
errU = sqrt(sum(wqJ(:).*err(:)))/sqrt(sum(wqJ(:).*uq(:)));

serr = (Vq*U{3} - sxx(xq,yq,FinalTime)).^2 + (Vq*U{4} - syy(xq,yq,FinalTime)).^2 + (Vq*U{5} - sxy(xq,yq,FinalTime)).^2;
sq = sxx(xq,yq,FinalTime).^2 + syy(xq,yq,FinalTime).^2 + sxy(xq,yq,FinalTime).^2;
errS = sqrt(sum(wqJ(:).*serr(:)))/sqrt(sum(wqJ(:).*sq(:)));

figure
pterr = abs(Vp*U{1} - v1(xp,yp,FinalTime)) + abs(Vp*U{2} - v2(xp,yp,FinalTime));
color_line3(xp,yp,pterr,pterr,'.')
title(sprintf('velocity L2 err = %g\n',errU))

figure
pterr = abs(Vp*U{3} - sxx(xp,yp,FinalTime)) + abs(Vp*U{4} - syy(xp,yp,FinalTime)) + abs(Vp*U{5} - sxy(xp,yp,FinalTime));
color_line3(xp,yp,pterr,pterr,'.')
title(sprintf('sigma L2 err = %g\n',errS))

L2err = sqrt(sum(wqJ(:).*(err(:)+serr(:))))/sqrt(sum(wqJ(:).*(uq(:)+sq(:))))
errU
errS



function [rhs] = ElasRHS2D(U,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

global Nfld mu lambda Vq Pq tau useWADG
global mapBd vmapBd mapBt vmapBt
global v1 v2

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

opt = 3;
if opt==1 % traction BCs
    nSx(mapB) = -2*(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -2*(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
elseif opt==2 % basic ABCs
    nSx(mapB) = -(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
    dU{1}(mapB) = -U{1}(vmapB);
    dU{2}(mapB) = -U{2}(vmapB);
elseif opt==3 % mixed BCs - traction left, dirichlet right
    nSx(mapBt) = -2*(nx(mapBt).*U{3}(vmapBt) + ny(mapBt).*U{5}(vmapBt));
    nSy(mapBt) = -2*(nx(mapBt).*U{5}(vmapBt) + ny(mapBt).*U{4}(vmapBt));
    dU{1}(mapBd) = 2*(v1(x(vmapBd),y(vmapBd),time)-U{1}(vmapBd));
    dU{2}(mapBd) = 2*(v2(x(vmapBd),y(vmapBd),time)-U{2}(vmapBd));
else
    
end

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
    flux{fld}(:) = fc{fld}(:) + tau{fld}.*fp{fld}(:) ;
end

% compute right hand sides of the PDE's
% rr{1} =  -divSx   +  LIFT*(Fscale.*flux{1})/2.0;
% rr{2} =  -divSy   +  LIFT*(Fscale.*flux{2})/2.0;
% rr{3} =  -du1dx   +  LIFT*(Fscale.*flux{3})/2.0;
% rr{4} =  -du2dy   +  LIFT*(Fscale.*flux{4})/2.0;
% rr{5} =  -du12dxy +  LIFT*(Fscale.*flux{5})/2.0;

rr{1} =  divSx   +  LIFT*(Fscale.*flux{1})/2.0;
rr{2} =  divSy   +  LIFT*(Fscale.*flux{2})/2.0;
rr{3} =  du1dx   +  LIFT*(Fscale.*flux{3})/2.0;
rr{4} =  du2dy   +  LIFT*(Fscale.*flux{4})/2.0;
rr{5} =  du12dxy +  LIFT*(Fscale.*flux{5})/2.0;


% C = [2*mu+lambda       lambda       0
%      lambda       2*mu+lambda       0
%           0       0                mu]

rhs{1} = rr{1};
rhs{2} = rr{2};
rhs{3} = (2*mu+lambda).*rr{3} + lambda.*rr{4};
rhs{4} = lambda.*rr{3} + (2*mu+lambda).*rr{4};
rhs{5} = (mu) .* rr{5};


return;

