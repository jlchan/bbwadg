function L2err = Elasticity2D_Stonely(N,K1D)

% clear all, clear
clear -global *

Globals2D

if nargin==0
    K1D = 4;
    N = 3;
end
c_flag = 0;
FinalTime = 1;

% filename = 'Grid/Other/block2.neu';
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
aa = 5;
% [Nv, VX, VY, K, EToV] = unif_tri_mesh(aa*K1D,K1D);
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D,aa*K1D);
VY = aa*VY;

StartUp2D;

% BuildPeriodicMaps2D(2,0);

% PlotMesh2D; axis on; return

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

[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
rp = .99*rp; sp = .99*sp;
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
% bfaces = sum(abs(abs(xfb)),1)<1e-8;
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
    hold on
    plot(x(vmapB),y(vmapB),'bo')
%     plot(x(vmapBt),y(vmapBt),'bo')
    
%     plot(x(vmapBd),y(vmapBd),'rd')
    plot(x,y,'k.')
    return
end
%%
global Nfld mu lambda rho Vq Pq tau useWADG
Nfld = 5; %(u1,u2,sxx,syy,sxy)

lambda1 = 3; mu1 = 3; rho1 = 10;
lambda2 = 1; mu2 = 1; rho2 = 1;

rho = rho2*ones(size(x));
lambda = lambda2*ones(size(x));
mu = mu2*ones(size(x));

useWADG = 1;
if useWADG
    rho = Vq*rho;
    mu = Vq*mu;
    lambda = Vq*lambda;
end

ids = yq > 0;
rho(ids) = rho1;
mu(ids) = mu1;
lambda(ids) = lambda1;
% plot3(x,y,Pq*rho,'.');return

% % piecewise constant
% mu = V\(Pq*mu); lambda = V\(Pq*lambda);
% lambda = repmat(lambda(1,:)*V(1,1),length(wq),1);
% mu = repmat(mu(1,:)*V(1,1),length(wq),1);
% % plot3(xq,yq,mu,'.');return

rho = Pq*rho;
mu = Pq*mu;
lambda = Pq*lambda;

tau0 = 1;
rhoF = max(rho(vmapM),rho(vmapP));
muF = max(mu(vmapM),mu(vmapP));
lambdaF = max(lambda(vmapM),lambda(vmapP));

rhoF = reshape(rhoF,Nfp*Nfaces,K);
muF = reshape(muF,Nfp*Nfaces,K);
lambdaF = reshape(lambdaF,Nfp*Nfaces,K);
tau{1} = tau0./rhoF;
tau{2} = tau0./rhoF;
tau{3} = tau0./(2*muF + lambdaF);
tau{4} = tau0./(2*muF + lambdaF);
tau{5} = tau0./(2*muF + lambdaF);
% for fld = 1:Nfld
%     tau{fld} = tau0*ones(size(muF));
% end
    

rho = Vq*rho;
mu = Vq*mu;
lambda = Vq*lambda;

%% params setup

% cp = sqrt((2*mu + lambda)/rho);
% cs = sqrt(mu/rho);

% assume constant
x0 = 0; y0 = .1;
p = exp(-16^2*((x-x0).^2 + (y-y0).^2));
u = zeros(Np,K);

U{1} = p;
U{2} = u;
U{3} = u;
U{4} = u;
U{5} = u;

cst = 0.546981324213884;
w = cst;
b1p = sqrt(1 - cst^2/((2*mu1+lambda1)/rho1));
b2p = sqrt(1 - cst^2/((2*mu2+lambda2)/rho2));
b1s = sqrt(1 - cst^2/(mu1/rho1));
b2s = sqrt(1 - cst^2/(mu2/rho2));
B1 = -1i * .2952173626624;
B2 = -.6798795208473;
B3 = 1i * 0.5220044931212;
B4 = -.9339639688697;

k = w/cst;

v1a = @(x,y,t) -w.*exp(k.*x.*1i - t.*w.*1i).*(B1.*k.*exp(-b1p.*k.*y).*1i + B2.*b1s.*k.*exp(-b1s.*k.*y)).*1i;
v2a = @(x,y,t) -w.*exp(k.*x.*1i - t.*w.*1i).*(B2.*k.*exp(-b1s.*k.*y).*1i - B1.*b1p.*k.*exp(-b1p.*k.*y)).*1i;
v1b = @(x,y,t) -k.*w.*exp(k.*x.*1i - t.*w.*1i).*(B3.*exp(b2p.*k.*y).*1i - B4.*b2s.*exp(b2s.*k.*y)).*1i;
v2b = @(x,y,t) -k.*w.*exp(k.*x.*1i - t.*w.*1i).*(B3.*b2p.*exp(b2p.*k.*y) + B4.*exp(b2s.*k.*y).*1i).*1i;

global v1 v2
v1 = @(x,y,t) real(v1a(x,y,t).*(y >= 0) + v1b(x,y,t).*(y < 0));
v2 = @(x,y,t) real(v2a(x,y,t).*(y >= 0) + v2b(x,y,t).*(y < 0));

u1ax   = @(x,y,t) k.*exp(k.*x.*1i - t.*w.*1i).*(B1.*k.*exp(-b1p.*k.*y).*1i + B2.*b1s.*k.*exp(-b1s.*k.*y)).*1i;
u2ay   = @(x,y,t) exp(k.*x.*1i - t.*w.*1i).*(B1.*b1p^2.*k^2.*exp(-b1p.*k.*y) - B2.*b1s.*k^2.*exp(-b1s.*k.*y).*1i);
u12axy = @(x,y,t) -exp(k.*x.*1i - t.*w.*1i).*(B2.*b1s^2.*k^2.*exp(-b1s.*k.*y) + B1.*b1p.*k^2.*exp(-b1p.*k.*y).*1i) ...
    + k.*exp(k.*x.*1i - t.*w.*1i).*(B2.*k.*exp(-b1s.*k.*y).*1i - B1.*b1p.*k.*exp(-b1p.*k.*y)).*1i;
u1bx   = @(x,y,t) k^2.*exp(k.*x.*1i - t.*w.*1i).*(B3.*exp(b2p.*k.*y).*1i - B4.*b2s.*exp(b2s.*k.*y)).*1i;
u2by   = @(x,y,t) k^2.*exp(k.*x.*1i - t.*w.*1i).*(B3.*b2p^2.*exp(b2p.*k.*y) + B4.*b2s.*exp(b2s.*k.*y).*1i);
u12bxy = @(x,y,t) -exp(k.*x.*1i - t.*w.*1i).*(B4.*b2s^2.*k^2.*exp(b2s.*k.*y) - B3.*b2p.*k^2.*exp(b2p.*k.*y).*1i) ...
    + k.*exp(k.*x.*1i - t.*w.*1i).*(B3.*b2p.*k.*exp(b2p.*k.*y) + B4.*k.*exp(b2s.*k.*y).*1i).*1i;

u1x = @(x,y,t) real(u1ax(x,y,t).*(y > 0) + u1bx(x,y,t).*(y < 0));
u2y = @(x,y,t) real(u2ay(x,y,t).*(y > 0) + u2by(x,y,t).*(y < 0));
u12xy = @(x,y,t) real(u12axy(x,y,t).*(y > 0) + u12bxy(x,y,t).*(y < 0));

fmu = @(x,y) mu1 * (y >= 0) + mu2*(y < 0);
flambda = @(x,y) lambda1 * (y >= 0) + lambda2*(y < 0);
sxx = @(x,y,t) ((2*fmu(x,y)+flambda(x,y)) .* u1x(x,y,t)) + (flambda(x,y).*u2y(x,y,t));
syy = @(x,y,t) flambda(x,y).*u1x(x,y,t) + ((2*fmu(x,y)+flambda(x,y)) .* u2y(x,y,t));
sxy = @(x,y,t) fmu(x,y) .* u12xy(x,y,t);

U{1} = Pq*v1(xq,yq,0);
U{2} = Pq*v2(xq,yq,0);
U{3} = Pq*((2*mu+lambda) .* u1x(xq,yq,0)) + Pq*(lambda.*u2y(xq,yq,0));
U{4} = Pq*(lambda.*u1x(xq,yq,0)) + Pq*((2*mu+lambda) .* u2y(xq,yq,0));
U{5} = Pq*(mu .* u12xy(xq,yq,0));

% for i = 3:5
%     U{i} = zeros(Np,K);
% end
% vv = Vp*(U{4}); color_line3(xp,yp,vv,vv,'.'); return

if 0
    [xp yp] = meshgrid(linspace(-4,4,200),linspace(-1,1,200));
    
    u1a = @(x,y,t) (1i.*k.*B1.*exp(-k.*b1p.*y) + k.*b1s.*B2.*exp(-k.*b1s.*y)).*exp(1i.*(k.*x-w.*t));
    u2a = @(x,y,t) (-k.*b1p.*B1.*exp(-k.*b1p.*y) + 1i.*k.*B2.*exp(-k.*b1s.*y)).*exp(1i.*(k.*x-w.*t));
    u1b = @(x,y,t) (1i.*k.*B3.*exp(k.*b2p.*y) - k.*b2s.*B4.*exp(k.*b2s.*y)).*exp(1i.*(k.*x-w.*t));
    u2b = @(x,y,t) (k.*b2p.*B3.*exp(k.*b2p.*y) + 1i.*k.*B4.*exp(k.*b2s.*y)).*exp(1i.*(k.*x-w.*t));
    
    for t = 0:.25:25
        sc = 2;
        clf
        ids = (yp > 0); xp1 = xp(ids); yp1 = yp(ids);
        plot(xp1-u1a(xp1,yp1,t)/sc,yp1+.01-u2a(xp1,yp1,t)/sc,'.')
        
        hold on
        
        ids = (yp < 0); xp2 = xp(ids); yp2 = yp(ids);
        plot(xp2-u1b(xp2,yp2,t)/sc,yp2-.01-u2b(xp2,yp2,t)/sc,'.')
        axis equal
%         axis([-2 2 -2 2 -1 1])
%         view(3)
        drawnow
        
    end
    return
end

%%

if 0
    tau0 = 1;
    
    rhoF = max(rho(vmapM),rho(vmapP));
    muF = max(mu(vmapM),mu(vmapP));
    lambdaF = max(lambda(vmapM),lambda(vmapP));
    
    rhoF = reshape(rhoF,Nfp*Nfaces,K);
    muF = reshape(muF,Nfp*Nfaces,K);
    lambdaF = reshape(lambdaF,Nfp*Nfaces,K);
    for fld = 1:5
        tau{fld} = tau0./rhoF;
        if fld > 2
            tau{fld} = tau0./(2*muF+lambdaF);
        end
    end
    
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
    plot(lam,'.','markersize',32)
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
CN = (N+1)^2/2; % guessing...
dt = .5/(sqrt(max((2*mu(:)+lambda(:))./rho(:)))*CN*max(Fscale(:)));

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
        p = U{3}+U{4}; % trace(S)
%         p = U{1};
        vp1 = Vp*p;        
        color_line3(xp,yp,vp1,vp1,'.');
        axis equal
        axis tight
        view(2)       
        
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

wqJ = diag(wq)*(Vq*J);
err = (Vq*U{1} - v1(xq,yq,FinalTime)).^2 + (Vq*U{2} - v2(xq,yq,FinalTime)).^2;
uq = v1(xq,yq,FinalTime).^2 + v2(xq,yq,FinalTime).^2;
errU = sqrt(sum(wqJ(:).*err(:)))/sqrt(sum(wqJ(:).*uq(:)));

serr = (Vq*U{3} - sxx(xq,yq,FinalTime)).^2 + (Vq*U{4} - syy(xq,yq,FinalTime)).^2 + (Vq*U{5} - sxy(xq,yq,FinalTime)).^2;
sq = sxx(xq,yq,FinalTime).^2 + syy(xq,yq,FinalTime).^2 + sxy(xq,yq,FinalTime).^2;
errS = sqrt(sum(wqJ(:).*serr(:)))/sqrt(sum(wqJ(:).*sq(:)));
if nargin==0
figure
pterr = abs(Vp*U{1} - v1(xp,yp,FinalTime)) + abs(Vp*U{2} - v2(xp,yp,FinalTime));
color_line3(xp,yp,pterr,pterr,'.')
title(sprintf('velocity L2 err = %g\n',errU))

figure
pterr = abs(Vp*U{3} - sxx(xp,yp,FinalTime)) + abs(Vp*U{4} - syy(xp,yp,FinalTime)) + abs(Vp*U{5} - sxy(xp,yp,FinalTime));
color_line3(xp,yp,pterr,pterr,'.')
title(sprintf('sigma L2 err = %g\n',errS))
end
L2err = sqrt(sum(wqJ(:).*(err(:)+serr(:))))/sqrt(sum(wqJ(:).*(uq(:)+sq(:))))
errU
errS



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

global v1 v2
opt = 3;
if opt==1 % traction BCs
    nSx(mapB) = -2*(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -2*(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
    
elseif opt==2 % basic ABCs
    nSx(mapB) = -(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
    dU{1}(mapB) = -U{1}(vmapB);
    dU{2}(mapB) = -U{2}(vmapB);
    
elseif opt==3
    
    dU{1}(mapB) = 2*(v1(x(vmapB),y(vmapB),time)-U{1}(vmapB));
    dU{2}(mapB) = 2*(v2(x(vmapB),y(vmapB),time)-U{2}(vmapB));
    
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
    flux{fld}(:) = fc{fld}(:) + tau{fld}(mapM).*fp{fld}(:);
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
    for fld = 1:Nfld
        rr{fld} = Vq*rr{fld};
    end
    rhs{1} = Pq*(rr{1}./rho);
    rhs{2} = Pq*(rr{2}./rho);
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

