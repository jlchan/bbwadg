function ElasticityInclusion2D

% clear all, clear
clear -global *

Globals2D

K1D = 8;
N = 4;
c_flag = 0;
FinalTime = .5;

% filename = 'Grid/Other/block2.neu';
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
VX = (VX+1)/2; 
VY = (VY+1)/2;

StartUp2D;

% BCType = zeros(size(EToV));
% for e = 1:K
%     BCType(e,EToE(e,:)==e) = 6;
% end
% 
% BuildBCMaps2D;
% nref = 1;
% for ref = 1:nref
%     Refine2D(ones(size(EToV)));
%     StartUp2D;
%     BuildBCMaps2D;
% end


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
global Nfld rho mu lambda Vq Pq tau useWADG
Nfld = 5; %(u1,u2,sxx,syy,sxy)

mu = ones(size(x));
lambda = 2*ones(size(x));
rho = ones(size(x));

mu0 = mu(1);
lambda0 = lambda(1);

% C = [2*mu0+lambda0       lambda0       0
%      lambda0       2*mu0+lambda0       0
%           0       0                mu0/2];

% mu = 1 + .5*sin(2*pi*x).*sin(2*pi*y);
% lambda = 1 + .5*sin(2*pi*x).*sin(2*pi*y);
% mu = 1 + .9*sin(8*pi*x).*sin(8*pi*y);
% lambda = 1 + .9*sin(8*pi*x).*sin(8*pi*y);

% ids =  mean(y) > .375 & mean(y) < .625 & mean(x) > .375 & mean(x) < .75;
% a = 0;
% mu(:,ids) = 5*mu(:,ids) + a*randn(Np,nnz(ids));
% lambda(:,ids) = 5*lambda(:,ids) + a*randn(Np,nnz(ids));

C = [2*mu(1) + lambda(1) lambda(1) 0;
    lambda(1) 2*mu(1)+lambda(1) 0;
    0 0 mu(1)];
% keyboard

useWADG = 1;
if useWADG
    rho = Vq*rho;
    mu = Vq*mu;
    lambda = Vq*lambda;
    if any(mu < 1e-8 | lambda < 1e-8)
        keyboard
    end
end

ids = yq > .375  & yq < .625 & xq > .375 & xq < .75;
mu(ids) = 5*mu(ids);
lambda(ids) = 5*lambda(ids);

tau0 = 1;
% for fld = 1:5
%     tau{fld} = tau0*ones(size(x));
%     if fld > 2
% %         tau{fld} = 1.5*tau0./max((2*mu(:)+lambda(:))./rho(:));
% %         tau{fld} = tau0./sqrt((2*mu(:)+lambda(:))./rho(:));
%         tau{fld} = tau0*ones(size(x));
%     end
% end
rhoF = max(rho(vmapM),rho(vmapP));
muF = max(mu(vmapM),mu(vmapP));
lambdaF = max(lambda(vmapM),lambda(vmapP));

rhoF = reshape(rhoF,Nfp*Nfaces,K);
muF = reshape(muF,Nfp*Nfaces,K);
lambdaF = reshape(lambdaF,Nfp*Nfaces,K);
tau{1} = tau0./rhoF;
tau{2} = tau0./rhoF;
tau{3} = 1.5*tau0./(2*muF + lambdaF);
tau{4} = 1.5*tau0./(2*muF + lambdaF);
tau{5} = 1.5*tau0./(2*muF + lambdaF);

%%
global mapBx vmapBx mapBt vmapBt t0 mapBL vmapBL

% t0 = .05;
t0 = .025;

% find x = 0 faces 
fbids = reshape(find(vmapP==vmapM),Nfp,nnz(vmapP==vmapM)/Nfp);
xfb = x(vmapP(fbids));
bfaces = sum(abs(xfb),1)<1e-8;
fbids = fbids(:,bfaces); fbids = fbids(:);
mapBL = fbids;
vmapBL = vmapP(mapBL);

fbids = reshape(find(vmapP==vmapM),Nfp,nnz(vmapP==vmapM)/Nfp);
xfb = x(vmapP(fbids));
bfaces = sum(abs(xfb),1)<1e-8 | sum(abs(xfb-1),1)<1e-8;
fbids = fbids(:,bfaces); fbids = fbids(:);
mapBx = fbids;
vmapBx = vmapP(mapBx);

% remove Dirichlet BCs for traction
mapBt = setdiff(mapB,mapBx);
vmapBt = vmapP(mapBt);

if 0
    plot(x,y,'.')
    hold on
    plot(x(vmapBx),y(vmapBx),'bo')
    plot(x(vmapBt),y(vmapBt),'rd')
    return
end

%% params setup

x0 = mean(VX); y0 = mean(VY) + .125;
x0 = .5; y0 = .5;
p = exp(-50^2*((x-x0).^2 + (y-y0).^2));
u = zeros(Np, K);

x0 = .1;
p = exp(-25^2*(x-x0).^2); % bar

% t0 = 0;
% U{1} = p; 
U{1} = u;
U{2} = u;
U{3} = u;
U{4} = u;
U{5} = u;

%% eig check

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
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
dt = 2/(max(2*mu(:)+lambda(:))*CN*max(Fscale(:)));
% CNh = (max(mu(:)+lambda(:))*CN*max(Fscale(:)));

% outer time step loop
tstep = 0;


figure
% colormap(gray)
colormap(hot)
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
        
        p = U{3}+U{4}; % trace(S)        
%         p = U{5};

        vv = Vp*p;
        vv = abs(vv);
        color_line3(xp,yp,vv,vv,'.');
        axis tight        
        caxis([.0,3])
        colorbar
        title(sprintf('time = %f',time)) 
        drawnow
%         view(3);        pause
        
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
    if mod(tstep,100)==0
        disp(sprintf('On timestep %d out of %d\n',tstep,round(FinalTime/dt)))
    end
end


p = U{3}+U{4}; % trace(S)
vv = abs(Vp*p);
color_line3(xp,yp,vv,vv,'.');
axis tight
title(sprintf('time = %f',time))


function [rhs] = ElasRHS2D(U,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

global Nfld rho mu lambda Vq Pq tau useWADG

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
 
% Impose left boundary conditions on velocity
global mapBx vmapBx mapBt vmapBt t0 mapBL vmapBL
if time < t0
    dU{1}(mapBL) = 2*(sin(pi*time/t0)*(time<t0) - U{1}(vmapBL));
    dU{2}(mapBL) = -2*U{2}(vmapBL);
else
    opt=3;
    if opt==1 % traction BCs        
        nSx(mapB) = -2*(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
        nSy(mapB) = -2*(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
        %     end
    elseif opt==2 % basic ABCs
        nSx(mapB) = -(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
        nSy(mapB) = -(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
        dU{1}(mapB) = -U{1}(vmapB);
        dU{2}(mapB) = -U{2}(vmapB);
    elseif opt==3 % mixed ABCs and traction
       
        nSx(mapBt) = -2*(nx(mapBt).*U{3}(vmapBt) + ny(mapBt).*U{5}(vmapBt));
        nSy(mapBt) = -2*(nx(mapBt).*U{5}(vmapBt) + ny(mapBt).*U{4}(vmapBt));
        
        nSx(mapBx) = -(nx(mapBx).*U{3}(vmapBx) + ny(mapBx).*U{5}(vmapBx));
        nSy(mapBx) = -(nx(mapBx).*U{5}(vmapBx) + ny(mapBx).*U{4}(vmapBx));
        dU{1}(mapBx) = -U{1}(vmapBx);
        dU{2}(mapBx) = -U{2}(vmapBx);

    end
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
%           0       0                mu/2]
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
    rhs{1} = rr{1}./rho;
    rhs{2} = rr{2}./rho;
    rhs{3} = (2*mu+lambda).*rr{3} + lambda.*rr{4};
    rhs{4} = lambda.*rr{3} + (2*mu+lambda).*rr{4};
    rhs{5} = (mu) .* rr{5};
end

return;


