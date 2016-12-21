function L2err = Elasticity2D_Lamb(N,K1D)

% clear all, clear
clear -global *

Globals2D

if nargin==0
    K1D = 16;
    N = 5;
end
c_flag = 0;
FinalTime = 1;

% filename = 'Grid/Other/block2.neu';
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
[Nv, VX, VY, K, EToV] = unif_tri_mesh(2*K1D,K1D);
VY = (VY+1)/2-.5;
% VX = 2*VX;

StartUp2D;

BuildPeriodicMaps2D(2,0);

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
% xfb = x(vmapM(fbids));
yfb = y(vmapM(fbids));
bfaces = sum(abs(abs(yfb)-.5),1)<1e-8;

% yfb = y(vmapM(fbids));
% bfaces = sum(abs(ptsfb),1)<1e-8;

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
global Nfld mu lambda rho Vq Pq tau useWADG
Nfld = 5; %(u1,u2,sxx,syy,sxy)

rho = 1;
lambda = 2;
mu = 1;

% global invC
% C = [2*mu+lambda       lambda       0
%      lambda       2*mu+lambda       0
%           0       0                mu];
% invC = inv(C);

tau0 = 1;
for fld = 1:5
    tau{fld} = tau0;
    if fld > 2
        %tau{fld} = tau0./(mu+lambda);
        tau{fld} = tau0*ones(size(mu));
    end
end

useWADG = 0;
if useWADG
    mu = Vq*mu;
    lambda = Vq*lambda;
end

%% check eigs


if 0
    tauvec = 0;
    for ii = 1:length(tauvec)
        for fld = 1:5
            tau{fld} = tauvec(ii);
            if fld > 2
                tau{fld} = tauvec(ii)./(mu+lambda);
            end
        end
        %         tau = tauvec(ii)*[1 1 1 1 1];
        u = zeros(Nfld*Np*K,1);
        rhs = zeros(Nfld*Np*K,1);
        A = zeros(Nfld*Np*K);
        ids = 1:Np*K;
        for i = 1:Nfld*Np*K
            u(i) = 1;
            for fld = 1:Nfld
                U{fld} = reshape(u(ids + (fld-1)*Np*K),Np,K);
            end
            rU = ElasRHS2D(U,0);
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
        points{ii} = [real(lam) imag(lam)];
        plot(lam,'.','markersize',32)
        hold on
        title(sprintf('Largest real part = %g\n',max(real(lam))))
        axis equal
        %         drawnow
        max(abs(lam))
    end
    keyboard
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


%%

time = 0;

% Runge-Kutta residual storage
for fld = 1:Nfld
    res{fld} = zeros(Np,K);
end

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
dt = 1/(max(mu(:)+lambda(:))*CN*max(Fscale(:)));

% outer time step loop
tstep = 0;


% figure
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
    
    if nargin==0 && mod(tstep,10)==0
        clf
        
        %         p = U{3}+U{4}; % trace(S)
        vp1 = Vp*U{1};
        vp1 = vp1/2.5e4;
        %         subplot(2,1,1)
        color_line3(xp,yp,vp1,vp1,'.');
        %         axis tight
        %         axis equal
        %         colorbar
        axis([-1 1 -.5 .5 -1.5 1.5])
        
        %         subplot(2,1,2)
        %         vv = u2(xp,yp,time)/2.5e3;
        % %         p = U{2};
        % %         vv = Vp*p;
        % %         vv = vv/1e3;
        %         color_line3(xp,yp,vv,vv,'.');
        % %         axis tight
        
        title(sprintf('time = %f',time))
        colorbar
        %         axis([-1 1 -.5 .5 -1.5 1.5])
        % axis equal
        %                 view(3);
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

L2err = sqrt(sum(wqJ(:).*(err(:)+serr(:))))/sqrt(sum(wqJ(:).*(uq(:)+sq(:))))
%  title(sprintf('L2 err = %g\n',L2err))

% vv = Vp*p;
% err = (abs(u1(xp,yp,FinalTime)-vv));
% figure
% color_line3(xp,yp,err,err,'.')




function [rhs] = ElasRHS2D(U,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

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

opt = 1;
if opt==1 % traction BCs
    nSx(mapB) = -2*(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -2*(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
    
elseif opt==2 % basic ABCs
    nSx(mapB) = -(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
    dU{1}(mapB) = -U{1}(vmapB);
    dU{2}(mapB) = -U{2}(vmapB);
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

return;

