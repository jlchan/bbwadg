function L2err = Elasticity2D(N,K1D)

% clear all, clear
clear -global *

Globals2D

if nargin==0
    K1D = 2;
    N = 5;
end
c_flag = 0;
FinalTime = 1;

% filename = 'Grid/Other/block2.neu';
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
VX = (VX+1)/2; 
VY = (VY+1)/2;
% VX = 2*VX;

StartUp2D;

BCType = zeros(size(EToV));
for e = 1:K
    BCType(e,EToE(e,:)==e) = 6;
end

BuildBCMaps2D;
nref = 1;
for ref = 1:nref
    Refine2D(ones(size(EToV)));
    StartUp2D;
    BuildBCMaps2D;
end


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

%%
global Nfld mu lambda Vq Pq tau useWADG
Nfld = 5; %(u1,u2,sxx,syy,sxy)

mu = 1;
lambda = 1;

mu0 = mu(1);
lambda0 = lambda(1);

% C = [2*mu0+lambda0       lambda0       0
%      lambda0       2*mu0+lambda0       0
%           0       0                mu0/2];

% mu = 1 + .5*sin(2*pi*x).*sin(2*pi*y);
% lambda = 1 + .5*sin(2*pi*x).*sin(2*pi*y);
% mu = 1 + .9*sin(8*pi*x).*sin(8*pi*y);
% lambda = 1 + .9*sin(8*pi*x).*sin(8*pi*y);

% ids =  mean(y) > .375 & mean(y) < .625 & mean(x) > .5 & mean(x) < .75;
% mu(:,ids) = 10*mu(:,ids);
% lambda(:,ids) = 10*lambda(:,ids);

% plot3(x,y,mu,'.');return
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


%% eigs

if 0
    tauvec = 0;
    for ii = 1:length(tauvec)
        for fld = 1:5
            tau{fld} = tauvec(ii)*ones(size(x));
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
            rU = ElasRHS2D(U);
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
k = 1; % frequency of solution
W = (2*k-1)/2*pi;

x0 = mean(VX); y0 = mean(VY) + .125;
x0 = .1; y0 = .5;
p = exp(-50^2*((x-x0).^2 + (y-y0).^2));
u = zeros(Np, K);

w = sqrt(2*mu);
u1 = @(x,y,t) cos(w*pi*t).*cos(pi*x).*sin(pi*y);
u2 = @(x,y,t) -cos(w*pi*t).*sin(pi*x).*cos(pi*y);

v1 = @(x,y,t) w.*pi.*sin(pi.*t.*w).*cos(pi.*x).*sin(pi.*y);
v2 = @(x,y,t) -w.*pi.*sin(pi.*t.*w).*cos(pi.*y).*sin(pi.*x);

u1x = @(x,y,t) -pi.*sin(pi.*x).*sin(pi.*y).*cos(pi.*t.*w);
u2y = @(x,y,t) pi.*sin(pi.*x).*sin(pi.*y).*cos(pi.*t.*w);
u12xy = @(x,y,t) zeros(size(x));

sxx = @(x,y,t) ((2*mu+lambda) .* u1x(x,y,t)) + (lambda.*u2y(x,y,t));
syy = @(x,y,t) lambda.*u1x(x,y,t) + ((2*mu+lambda) .* u2y(x,y,t));
sxy = @(x,y,t) mu .* u12xy(x,y,t);

U{1} = v1(x,y,0);
U{2} = v2(x,y,0);
U{3} = sxx(x,y,0);
U{4} = syy(x,y,0);
U{5} = sxy(x,y,0);


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
CNh = (max(mu(:)+lambda(:))*CN*max(Fscale(:)));

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
        
        p = U{5}; % trace(S)        
%         p = U{1};

        vv = Vp*p;        
%         vv = abs(vv);
%         vv = max(vv(:))-vv;
        color_line3(xp,yp,vv,vv,'.');
        axis tight        
        title(sprintf('time = %f',time)) 
        view(3);
        axis([0 1 0 1 -5 5])
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

% figure
% pterr = abs(Vp*U{1} - v1(xp,yp,FinalTime)) + abs(Vp*U{2} - v2(xp,yp,FinalTime));
% color_line3(xp,yp,pterr,pterr,'.')
% title(sprintf('velocity L2 err = %g\n',errU))
% 
% figure
% pterr = abs(Vp*U{3} - sxx(xp,yp,FinalTime)) + abs(Vp*U{4} - syy(xp,yp,FinalTime)) + abs(Vp*U{5} - sxy(xp,yp,FinalTime));
% color_line3(xp,yp,pterr,pterr,'.')
% title(sprintf('sigma L2 err = %g\n',errS))

L2err = sqrt(sum(wqJ(:).*(err(:)+serr(:))))/sqrt(sum(wqJ(:).*(uq(:)+sq(:))))
if nargin==0
    title(sprintf('L2 err = %g\n',L2err))
end
% vv = Vp*p;
% err = (abs(u1(xp,yp,FinalTime)-vv));
% figure
% color_line3(xp,yp,err,err,'.')



function [rhs] = ElasRHS2D(U,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

global Nfld mu lambda Vq Pq tau useWADG

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

nSx = nx.*dU{3} + ny.*dU{5};
nSy = nx.*dU{5} + ny.*dU{4};
nUx = dU{1}.*nx;
nUy = dU{2}.*ny;
nUxy = dU{2}.*nx + dU{1}.*ny;

% dU{1}(mapB) = -2*U{1}(vmapBx);
% dU{2}(mapB) = -2*U{2}(vmapB);
% dU{5}(mapB) = -2*U{5}(vmapB);

% evaluate upwind fluxes
fc{1} = nSx;
fc{2} = nSy;
fc{3} = nUx;
fc{4} = nUy;
fc{5} = nUxy;

% penalization terms
fp{1} = dU{1} + nx.*ny.*dU{2};
fp{2} = dU{2} + nx.*ny.*dU{1};
fp{3} = nx.*nSx;
fp{4} = ny.*nSy;
fp{5} = dU{5} + nx.*ny.*(dU{3} + dU{4});

% zero traction BCs - set nx*Sxx + ny*Sxy = nx*Sxy + ny*Syy = 0
nSxM = nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB);
nSyM = nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB);
fc{1}(mapB) = -2*nSxM;
fc{2}(mapB) = -2*nSyM;
fp{3}(mapB) = -2*nx(mapB).*nSxM;
fp{4}(mapB) = -2*ny(mapB).*nSyM;

% global invC
% tmp{1} = invC(1,1)*fp{3} + invC(1,2)*fp{4};
% tmp{2} = invC(2,1)*fp{3} + invC(2,2)*fp{4};
% tmp{3} = invC(3,3)*fp{5};
% fp{3} = tmp{1};
% fp{4} = tmp{2};
% fp{5} = tmp{3};

flux = cell(5,1);
for fld = 1:Nfld   
    flux{fld} = zeros(Nfp*Nfaces,K);
    flux{fld}(:) = tau{fld}*fp{fld}(:) - fc{fld}(:);
end


% compute right hand sides of the PDE's
rr{1} =  -divSx   +  LIFT*(Fscale.*flux{1})/2.0;
rr{2} =  -divSy   +  LIFT*(Fscale.*flux{2})/2.0;
rr{3} =  -du1dx   +  LIFT*(Fscale.*flux{3})/2.0;
rr{4} =  -du2dy   +  LIFT*(Fscale.*flux{4})/2.0;
rr{5} =  -du12dxy +  LIFT*(Fscale.*flux{5})/2.0;

% C = [2*mu+lambda       lambda       0
%      lambda       2*mu+lambda       0
%           0       0                mu/2]

if useWADG
    for fld = 3:Nfld
        rr{fld} = Vq*rr{fld};
    end
    rhs{1} = rr{1};
    rhs{2} = rr{2};
    rhs{3} = Pq*((2*mu+lambda).*rr{3} + lambda.*rr{4});
    rhs{4} = Pq*(lambda.*rr{3} + (2*mu+lambda).*rr{4});
    rhs{5} = Pq*((mu) .* rr{5});
else
    rhs{1} = rr{1};
    rhs{2} = rr{2};
    rhs{3} = (2*mu+lambda).*rr{3} + lambda.*rr{4};
    rhs{4} = lambda.*rr{3} + (2*mu+lambda).*rr{4};
    rhs{5} = (mu) .* rr{5};
end

return;

