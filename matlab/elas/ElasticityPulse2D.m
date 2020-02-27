function ElasticityPulse2D

% clear all, clear
clear -global *

Globals2D

K1D = 8;
N = 5;
c_flag = 0;
FinalTime = 2.0;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

StartUp2D;

[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

Nq = 2*N;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
Pq = V*V'*Vq'*diag(wq); % J's cancel out
Mref = inv(V*V');
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;


%%
global Nfld mu lambda Vq Pq tau useWADG
Nfld = 5; %(u1,u2,sxx,syy,sxy)

mu = ones(size(x));
lambda = ones(size(x));

mu0 = mu(1);
lambda0 = lambda(1);

% C = [2*mu0+lambda0       lambda0       0
%      lambda0       2*mu0+lambda0       0
%           0       0                mu0];

% mu = 1 + .5*sin(2*pi*x).*sin(2*pi*y);
% lambda = 1 + .5*sin(2*pi*x).*sin(2*pi*y);

% ids =  mean(y) > .375 & mean(y) < .625 & mean(x) > .5 & mean(x) < .75;
% ids = mean(y) > 0;
% mu(:,ids) = 0*mu(:,ids);
% lambda(:,ids) = lambda(:,ids);


global C11 C12 C13 C22 C23 C33

useWADG = 0;
if useWADG
    mu = Vq*mu;
    lambda = Vq*lambda;   
end

tau0 = 1;
for fld = 1:5
    tau{fld} = tau0*ones(size(x));
    if fld > 2
%         tau{fld} = tau0./(2*mu+lambda);
        tau{fld} = tau0*ones(size(x));
%         tau{fld} = tau0./Cnorm;
%         tau{fld} = tau0./max(Cnorm(:))*ones(size(x));
    end
end


%% params setup

x0 = 0;
y0 = 0;
p = exp(-10^2*((x-x0).^2 + (y-y0).^2));
% p = exp(-50^2*((x-x0).^2));
u = zeros(Np, K);

% p = (1-x).*(1+x).*y.*(1-y).*(y>0);
U{1} = u;
U{2} = p;
U{3} = u;
U{4} = u;
U{5} = u;

%%
if 0
    
    t0 = 0;
    
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
    
    if 1 && mod(tstep,10)==0
        clf
        
        p = (U{3} + U{4})/2; % trace(S)        
        vv = Vp*p;
        color_line3(xp,yp,vv,vv,'.');
        axis tight        
        title(sprintf('time = %f',time));
        colorbar;
        caxis([-.1,.1])

        drawnow
                
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
    
    if mod(tstep,100)==0
        disp(sprintf('On timestep %d out of %d\n',tstep,round(FinalTime/dt)))
    end
end

clf
p = (U{3} + U{4})/2; % trace(S)
vv = Vp*p;
color_line3(xp,yp,vv,vv,'.');
axis tight
title(sprintf('time = %f',time));
colorbar;

pelas = p;

keyboard



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
 
opt=2;
if opt==1 % traction BCs
    %     global mapBL vmapBL t0
    %     if time < t0
    %         dU{1}(mapBL) = 2*(sin(pi*time/t0)*(time<t0) - U{1}(vmapBL));
    %         dU{2}(mapBL) = -2*U{2}(vmapBL);
    %     else
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
%     rhs{3} = Pq*(C11.*rr{3} + C12.*rr{4} + C13.*rr{5});
%     rhs{4} = Pq*(C12.*rr{3} + C22.*rr{4} + C23.*rr{5});
%     rhs{5} = Pq*(C13.*rr{3} + C23.*rr{4} + C33.*rr{5});
else
    rhs{1} = rr{1};
    rhs{2} = rr{2};
    rhs{3} = (2*mu+lambda).*rr{3} + lambda.*rr{4};
    rhs{4} = lambda.*rr{3} + (2*mu+lambda).*rr{4};
    rhs{5} = (mu) .* rr{5};
end

return;

