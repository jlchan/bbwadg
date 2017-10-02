function rel_err = ElasticityReferenceSolution2D(N,K1D)

% clear all, clear
clear -global *

Globals2D
if nargin==0
K1D = 32;
N = 5;
end
c_flag = 0;
FinalTime = .5;

filename = 'Grid/Other/block2.neu';
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
% VX = (VX+1)/2; 
% VY = (VY+1)/2;
% VX = 2*VX;

StartUp2D;

% BuildPeriodicMaps2D(2,2)

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
   
    a = .25; b = 1;
    mu = 1 + a*cos(b*pi*xq).*cos(b*pi*yq);
    lambda = 1 + a*sin(b*pi*xq).*sin(b*pi*yq);
    
else
    mu = ones(size(x));
    lambda = ones(size(x));
    
end

tau0 = 1;
for fld = 1:5
    tau{fld} = tau0*ones(size(x));
    if fld > 2
%         tau{fld} = tau0.*ones(size(x))/max(max(2*(Pq*mu)+(Pq*lambda)));
        tau{fld} = tau0*ones(size(x));
    end
end


%% params setup

mu0 = 1;
w = sqrt(2*mu0);
u1 = @(x,y,t) cos(w*pi*t).*cos(pi*x).*sin(pi*y);
u2 = @(x,y,t) -cos(w*pi*t).*sin(pi*x).*cos(pi*y);

v1 = @(x,y,t) w.*pi.*sin(pi.*t.*w).*cos(pi.*x).*sin(pi.*y);
v2 = @(x,y,t) -w.*pi.*sin(pi.*t.*w).*cos(pi.*y).*sin(pi.*x);
% p = exp(-100*(x.^2+y.^2));
u = zeros(Np, K);

U{1} = u1(x,y,0);
U{2} = u2(x,y,0);
U{3} = u;
U{4} = u;
U{5} = u;


%%

time = 0;

% Runge-Kutta residual storage
for fld = 1:Nfld
    res{fld} = zeros(Np,K);
end

% compute time step size
CN = (N+1)^2/2; % guessing...
dt = 2/(max(2*mu(:)+lambda(:))*CN*max(Fscale(:)));

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
        
        p = U{3}; % trace(S)        
%         p = U{1};

        vv = Vp*p;
%         vv = abs(vv);
%         vv = max(vv(:))-vv;
        color_line3(xp,yp,vv,vv,'.');
        axis tight        
        title(sprintf('time = %f',time));
        colorbar;
        view(2); 
        axis equal
        axis([-1 1 -1 1 -2 2])
        drawnow
        
        
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
    
%     Mu = M*(J.*U{1});
%     uu = U{1}(:)'*Mu(:);
%     
%     Mu = M*(J.*U{2});
%     uu = uu + U{2}(:)'*Mu(:);
%     
%     u1q = Vq*(U{3});
%     u2q = Vq*(U{4});
%     u3q = Vq*(U{5});
%         
%     Cu1 = iC11.*u1q + iC12.*u2q + iC13.*u3q;
%     Cu2 = iC12.*u1q + iC22.*u2q + iC23.*u3q;
%     Cu3 = iC13.*u1q + iC23.*u2q + iC33.*u3q;
%     uu = uu + sum(sum(wqJ.*(u1q.*Cu1 + u2q.*Cu2 + u3q.*Cu3)));
%         
%     unorm2(tstep) = uu;
%     tvec(tstep) = time;
    
    if mod(tstep,100)==0
        disp(sprintf('On timestep %d out of %d\n',tstep,round(FinalTime/dt)))
    end
end

if nargin==0
    clf
    p = U{3};% + U{4}; % trace(S)
    vv = Vp*p;
    color_line3(xp,yp,vv,vv,'.');
    axis tight
    title(sprintf('time = %f',time));
    colorbar;
    view(3);
    axis([-1 1 -1 1 -2 2])
end

Nq = 3*N;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
Pq = V*V'*Vq'*diag(wq); % J's cancel out
Mref = inv(V*V');
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;
wqJ = diag(wq)*Jq;

Nref = 50;
[Uref, invVquad] = ElasticitySpectral2D(Nref,FinalTime);
Vqquad = Vandermonde2DQuad(Nref,xq(:),yq(:))*invVquad;

err = 0;
Unorm = 0;
for fld = 1:Nfld    
    Uq1 = Vq*U{fld};
    Uq2 = reshape(Vqquad * Uref{fld}(:),length(wq),K);
    err = err + sum(sum(wqJ.*(Uq1-Uq2).^2));
    Unorm = Unorm + sum(sum(wqJ.*(Uq2).^2));    
end

rel_err = sqrt(err) / sqrt(Unorm)

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
 
opt=1;
if opt==1 % traction BCs    
    nSx(mapB) = -2*(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -2*(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
    %     end    
elseif opt==2 % basic ABCs
    nSx(mapB) = -(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
    dU{1}(mapB) = -U{1}(vmapB);
    dU{2}(mapB) = -U{2}(vmapB);    
elseif opt==3 % zero velocity
    dU{1}(mapB) = -2*U{1}(vmapB);
    dU{2}(mapB) = -2*U{2}(vmapB);    
    
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

