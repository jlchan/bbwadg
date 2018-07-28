%% make symbolic fxns for manufactured soln
% syms k x y cP cS t a b
% 
% u1 = cos(k*((a*x+b*y)-cP*t));
% u2 = cos(k*((a*x+b*y)-cS*t));
% 
% v1 = matlabFunction(simplify(diff(u1,t)))
% v2 = matlabFunction(simplify(diff(u2,t)))
% 
% u1x = matlabFunction(simplify(diff(u1,x)))
% u2y = matlabFunction(simplify(diff(u2,y)))
% u12xy = matlabFunction(simplify(diff(u1,y) + diff(u2,x)))
% 
% 
% v1x = matlabFunction(diff(simplify(diff(u1,t)),x))
% v2y = matlabFunction(diff(simplify(diff(u2,t)),y))

function L2err = ElasticityManufacturedSolution2D(N,K1D,tau0)

% clear all, clear
clear -global *

Globals2D
if nargin==0
    K1D = 4;
    N = 5;
end
if nargin~=3
    tau0=1;
end

FinalTime = 5; %2*sqrt(2);

if 0
    if K1D==2
        filename = 'Grid/Maxwell2D/Maxwell1.neu';
    elseif K1D==4
        filename = 'Grid/Maxwell2D/Maxwell05.neu';
    elseif K1D==8
        filename = 'Grid/Maxwell2D/Maxwell025.neu';
    elseif K1D==16
        filename = 'Grid/Maxwell2D/Maxwell0125.neu';
    elseif K1D==32
        filename = 'Grid/Maxwell2D/Maxwell00625.neu';
    end
    [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
else
    [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
end

StartUp2D;

BuildPeriodicMaps2D(2,2)

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

useWADG = 1;
if useWADG
    
    mu = 1;
    
    a = 1; % variation in lambda
    w = 2;
    lambda0 = 2;
    lamOsc = @(x,y) a*.5*sin(w*pi*x).*sin(w*pi*y);
    lambdafun = @(x,y) lambda0 + lamOsc(x,y);    
    lambda = lambdafun(xq,yq); 

%     vv = lambdafun(xp,yp); color_line3(xp,yp,vv,vv,'.'); return
    
else
    mu = ones(size(x));
    lambda = ones(size(x));
    
end

for fld = 1:5
    tau{fld} = tau0*ones(size(x));    
end


%% params setup

beta = [ 1 0];
a = beta(1)/norm(beta); 
b = beta(2)/norm(beta);
k = pi;
cP = @(x,y) sqrt(2*mu + lambda0);
cS = @(x,y) sqrt(mu);


v1  = @(x,y,t) cP(x,y).*k.*sin(k.*(-cP(x,y).*t + a.*x + b.*y));
v2  = @(x,y,t) cS(x,y).*k.*sin(k.*(-cS(x,y).*t + a.*x + b.*y));
                    
u1x = @(x,y,t) -a.*k.*sin(k.*(-cP(x,y).*t+a.*x+b.*y));
u2y = @(x,y,t) -b.*k.*sin(k.*(-cS(x,y).*t+a.*x+b.*y));
u12xy = @(x,y,t) -a.*k.*sin(k.*(-cS(x,y).*t+a.*x+b.*y))-b.*k.*sin(k.*(-cP(x,y).*t+a.*x+b.*y));
sxx = @(x,y,t) ((2*mu+lambda0) .* u1x(x,y,t)) + (lambda0.*u2y(x,y,t));
syy = @(x,y,t) lambda0.*u1x(x,y,t) + ((2*mu+lambda0) .* u2y(x,y,t));
sxy = @(x,y,t) mu .* u12xy(x,y,t);
  
U{1} = Pq*(v1(xq,yq,0));
U{2} = Pq*(v2(xq,yq,0));
U{3} = Pq*(sxx(xq,yq,0));
U{4} = Pq*(syy(xq,yq,0));
U{5} = Pq*(sxy(xq,yq,0));

% source terms for varying wavespeed
global fsrc xq yq
v1x = @(x,y,t) cP(x,y).*k.^2.*cos(k.*(x-cP(x,y).*t));
v2y = @(x,y,t) 0*x;
fsrc = @(x,y,t) lamOsc(x,y).*(v1x(x,y,t)+v2y(x,y,t));

% for t = 0:.01:5
%     clf
%     vv = SS(xp,yp,t)+2*SP(xp,yp,t);
%     color_line3(xp,yp,vv,vv,'.')
%     caxis([-2,2])
%     
%     drawnow    
% end
% return

% keyboard

%%

time = 0;

% Runge-Kutta residual storage
for fld = 1:Nfld
    res{fld} = zeros(Np,K);
end

% compute time step size
CN = (N+1)^2/2; % guessing...
dt = 1/(max(2*mu(:)+lambda(:))*CN*max(Fscale(:)));

% outer time step loop
tstep = 0;

M = inv(V*V');
wqJ = diag(wq)*(Vq*J);

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
        
        subplot(1,2,1)
        p = U{2};
        vv = Vp*p; 
        color_line3(xp,yp,vv,vv,'.');
        axis tight
        title(sprintf('time = %f',time));
        colorbar;
        caxis([-10,10])
        view(2);
%         zlim([-3 3])
%         axis equal
        
        subplot(1,2,2)
%         p = U{4};
%         vv = Vp*p;
        vv = v2(xp,yp,time);
        color_line3(xp,yp,vv,vv,'.');        
        axis tight
        title(sprintf('time = %f',time));
        colorbar;
        caxis([-10,10])
        view(2);
%         zlim([-3 3])
%         axis equal
        
        drawnow               
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
        
    if mod(tstep,100)==0
        disp(sprintf('On timestep %d out of %d\n',tstep,round(FinalTime/dt)))
    end
end

% if nargin==0
%     clf
%     p = U{3};% + U{4}; % trace(S)
%     vv = Vp*p;
%     color_line3(xp,yp,vv,vv,'.');
%     axis tight
%     title(sprintf('time = %f',time));
%     colorbar;
%     view(2);
% end

Nq = 3*N;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
Pq = V*V'*Vq'*diag(wq); % J's cancel out
Mref = inv(V*V');
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;
wqJ = diag(wq)*Jq;

wqJ = diag(wq)*(Vq*J);
err = (Vq*U{1} - v1(xq,yq,FinalTime)).^2 + (Vq*U{2} - v2(xq,yq,FinalTime)).^2;
uq = v1(xq,yq,FinalTime).^2 + v2(xq,yq,FinalTime).^2;
errU = sqrt(sum(wqJ(:).*err(:)))/sqrt(sum(wqJ(:).*uq(:)));

serr = (Vq*U{3} - sxx(xq,yq,FinalTime)).^2 + (Vq*U{4} - syy(xq,yq,FinalTime)).^2 + (Vq*U{5} - sxy(xq,yq,FinalTime)).^2;
sq = sxx(xq,yq,FinalTime).^2 + syy(xq,yq,FinalTime).^2 + sxy(xq,yq,FinalTime).^2;
errS = sqrt(sum(wqJ(:).*serr(:)))/sqrt(sum(wqJ(:).*sq(:)));

L2err = sqrt(sum(wqJ(:).*(err(:)+serr(:))))/sqrt(sum(wqJ(:).*(uq(:)+sq(:))))

% keyboard



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

% velocity fluxes
nSx = nx.*dU{3} + ny.*dU{5};
nSy = nx.*dU{5} + ny.*dU{4};

% opt=1;
% if opt==1 % traction BCs
%     nSx(mapB) = -2*(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
%     nSy(mapB) = -2*(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
%     %     end
% elseif opt==2 % basic ABCs
%     nSx(mapB) = -(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
%     nSy(mapB) = -(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
%     dU{1}(mapB) = -U{1}(vmapB);
%     dU{2}(mapB) = -U{2}(vmapB);
% elseif opt==3 % zero velocity
%     dU{1}(mapB) = -2*U{1}(vmapB);
%     dU{2}(mapB) = -2*U{2}(vmapB);   
% end

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
    flux{fld}(:) = fc{fld}(:) + tau{fld}(vmapM).*fp{fld}(:);
end

% compute right hand sides of the PDE's
rr{1} =  divSx   +  LIFT*(Fscale.*flux{1})/2.0;
rr{2} =  divSy   +  LIFT*(Fscale.*flux{2})/2.0;
rr{3} =  du1dx   +  LIFT*(Fscale.*flux{3})/2.0;
rr{4} =  du2dy   +  LIFT*(Fscale.*flux{4})/2.0;
rr{5} =  du12dxy +  LIFT*(Fscale.*flux{5})/2.0;

% C = [2*mu+lambda       lambda       0
%      lambda       2*mu+lambda       0
%           0       0                mu/2]
for fld = 3:Nfld
    rr{fld} = Vq*rr{fld};
end
rhs{1} = rr{1};
rhs{2} = rr{2};

global fsrc xq yq
S = -Pq*fsrc(xq,yq,time);

rhs{3} = Pq*((2*mu+lambda).*rr{3} + lambda.*rr{4}) + S;
rhs{4} = Pq*(lambda.*rr{3} + (2*mu+lambda).*rr{4}) + S;
rhs{5} = Pq*((mu) .* rr{5});


return;

