% this code tests trace vs error based PAT iteration.  observations:
% - trace-based recon does better for about 1 iteration but error increases too fast
% - error-based recon doesn't increase error but doesn't improve results much
% - all results pretty sensitive to final time
%
% - try ben's code
% - add filter


function Wave2D

% clear all, clear
clear -global *

Globals2D

N = 3;
K1D = 8;
FinalTime = 2.1;

plotFlag = 0;
useTraces = 1;

global tau0
tau0 = 1;

cfun = @(x,y) ones(size(x));
% cfun = @(x,y) 1 + (x > 0);
% cfun = @(x,y) 1 + .25*sin(2*pi*x).*sin(2*pi*y); % smooth velocity
% cfun = @(x,y) (1 + .25*sin(2*pi*x).*sin(2*pi*y) + (y > 0)); % piecewise smooth velocity

% filename = 'Grid/Other/block2.neu';
% filename = 'Grid/Maxwell2D/Maxwell05.neu';
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
StartUp2D;

global Pq cq Vq

% PlotMesh2D; return

[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
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
cq = cfun(xq,yq);

% clf; vv = Vp*Pq*cq; color_line3(xp,yp,vv,vv,'.');axis equal;axis tight;return

% FinalTime = 2/min(cq(:));

%% params setup

x0 = .25; y0 = 0;

% gaussian pulse
pex = @(x,y) exp(-5^2*((x-x0).^2 + (y-y0).^2));
% pex = @(x,y) cos(.5*pi*x).*cos(.5*pi*y);

% % smooth annulus
% sig = @(x) 1-1./(1+exp(-100*x));
% r2 = @(x,y) x.^2 + y.^2;
% pex = @(x,y) sig((r2(x,y)-.25).^2);

% % boxes overlapping
% a = 2/K1D;
% pex = @(x,y) (abs(x+a)<.5 & abs(y)<.5) + (abs(x)<.5 & abs(y-a)<.5) + (abs(x-a)<.5 & abs(y+a)<.5);

% projection
p = Pq*pex(xq,yq);
u = zeros(Np, K);
v = zeros(Np, K);

% vv = Vp*p; color_line3(xp,yp,vv,vv,'.'); axis equal; axis tight; return

%% check eigs

if 0 & 3*Np*K < 4000
    U = zeros(Np*K,3); A = zeros(3*Np*K);
    for i = 1:3*Np*K
        U(i) = 1;
        p = reshape(U(:,1),Np,K); u = reshape(U(:,2),Np,K); v = reshape(U(:,3),Np,K);
        [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,0);
        A(:,i) = [rhsp(:); rhsu(:); rhsv(:)];
        U(i) = 0;
        if mod(i,round(3*Np*K/10))==0
            disp(sprintf('i = %d out of %d\n',i,3*Np*K))
        end
    end
    lam = eig(A);
    keyboard
end

%% setup timestepper

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
CNh = max(max(cq(:))*CN*max(Fscale(:)));
dt = 1.5/CNh;

Nstep = ceil(FinalTime/dt);
dt = FinalTime/Nstep;


%% generate synthetic data

NstepRK = Nstep*5;
pbc = zeros(length(mapB(:)),NstepRK);

% outer time step loop
for i = 1:Nstep
    for INTRK = 1:5
        [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,dt,'abc');
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        p = p+rk4b(INTRK)*resp;
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        
        % save traces
        id0 = ((i-1)*5+(INTRK-1)); % zero index
        pbc(:,NstepRK - id0) = p(vmapB);        
    end;
    
    if plotFlag && mod(i,10)==0
        vv = Vp*p;
        clf; color_line3(xp,yp,vv,vv,'.');
        axis equal; axis tight; colorbar
        title(sprintf('time = %f, max solution val = %f',i*dt,max(abs(vv(:)))))
        drawnow
    end
    if (mod(i,ceil(Nstep/10))==0)
        disp(sprintf('generating synthetic data: tstep %d out of %d\n',i,Nstep))
    end
end

disp('Synthetic data generated. ')


%% playback : reverse propagate to generate Rm from pressure boundary data

dt = -dt;
p = zeros(Np,K); u = p; v = p;
resp = zeros(Np,K); resu = resp; resv = resp;

for i = 1:Nstep
    for INTRK = 1:5
        id0 = ((i-1)*5+(INTRK-1)); % 1-index
        
        [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,dt,'rd',pbc(:,id0+1));
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        p = p + rk4b(INTRK)*resp;
        u = u + rk4b(INTRK)*resu;
        v = v + rk4b(INTRK)*resv;
    end;
    
    if plotFlag && mod(i,10)==0
        vv = Vp*p;
        clf; color_line3(xp,yp,vv,vv,'.');
        axis equal; axis tight; colorbar
        title(sprintf('time = %f, max solution val = %f',i*dt,max(abs(vv(:)))))
        drawnow
    end
    if (mod(i,ceil(Nstep/5))==0)
        disp(sprintf('first backwards prop: tstep %d out of %d',i,Nstep))
    end    
end

% save reversed data
p0 = p;

% plot initial cond
p_recon = p;
p0ex = Vq*Pq*pex(xq,yq); % use projection instead of exact?
err = diag(wq)*(Vq*J).*(Vq*p_recon-p0ex).^2;
L2err(1) = sqrt(sum(err(:)));
uq = diag(wq)*(Vq*J).*pex(xq,yq).^2;
Unorm = sqrt(sum(uq(:)));
disp(['Time reversal initial recon: L2 err = ',num2str(L2err/Unorm)])

%% iterate

% initialize to A*L*f
p = p0;

p_recon = p0;

for iter = 1:10
        
    p_old = p;
    
    % Assume p = (I - A*L)*p;    
    u = zeros(Np,K); v = zeros(Np,K);
    resp = zeros(Np,K); resu = resp; resv = resp;
    
    % apply L: fwd propagation and recover boundary measurements
    dt = abs(dt);
    for i = 1:Nstep
        for INTRK = 1:5            
            [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,dt,'abc');
            resp = rk4a(INTRK)*resp + dt*rhsp;
            resu = rk4a(INTRK)*resu + dt*rhsu;
            resv = rk4a(INTRK)*resv + dt*rhsv;
            p = p + rk4b(INTRK)*resp;
            u = u + rk4b(INTRK)*resu;
            v = v + rk4b(INTRK)*resv;
            
            if useTraces
                % save traces
                id0 = ((i-1)*5+(INTRK-1)); % zero index
                pbc(:,NstepRK - id0) = p(vmapB);
            end
        end;
        
        if plotFlag && mod(i,10)==0
            vv = Vp*p;
            clf; color_line3(xp,yp,vv,vv,'.');
            axis equal; axis tight; colorbar
            title(sprintf('time = %f, max solution val = %f',i*dt,max(abs(vv(:)))))
            drawnow
        end
        if (mod(i,ceil(Nstep/5))==0)
            disp(sprintf('fwd prop: tstep %d out of %d',i,Nstep))
        end
        
    end
    
    % compute new auxiliary initial conditions    
    if useTraces
        if 1 % compute harmonic extension: does way worse?
            [R vmapBT] = getCGRestriction();
            Mhat = inv(V*V');
            M = kron(spdiag(J(1,:)),Mhat);
            Dx = kron(spdiag(rx(1,:)),Dr) + kron(spdiag(sx(1,:)),Ds);
            Dy = kron(spdiag(ry(1,:)),Dr) + kron(spdiag(sy(1,:)),Ds);
            KK = Dx'*M*Dx + Dy'*M*Dy;
            KK = R*KK*R';
            KK(vmapBT) = speye(size(vmapBT));
            pB = zeros(Np,K); 
            pB(vmapB) = pbc(:,1);            
            pCG = diag(1./sum(R,2))*R*pB(:); % nodal averaging for BCs
            b = zeros(size(KK,2),1);
            b(vmapBT) = pCG(vmapBT);
            p = reshape(R'*(KK\b),Np,K);
%             clf
%             vv= Vp*p;
%             color_line3(xp,yp,vv,vv,'.')
%             return
        else
            p = zeros(size(p)); 
        end
        u = zeros(Np,K);
        v = zeros(Np,K);                    
    end
    
    % apply A: compute backwards propagation from boundary measurements, recover initial condition
    dt = -dt;
    for i = 1:Nstep
        for INTRK = 1:5
            if useTraces
                id0 = ((i-1)*5+(INTRK-1)); % 1-index
                [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,dt,'rd',pbc(:,id0+1));
            else
                [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,dt,'d');
            end
            
            resp = rk4a(INTRK)*resp + dt*rhsp;
            resu = rk4a(INTRK)*resu + dt*rhsu;
            resv = rk4a(INTRK)*resv + dt*rhsv;
            p = p + rk4b(INTRK)*resp;
            u = u + rk4b(INTRK)*resu;
            v = v + rk4b(INTRK)*resv;
        end;
        
        if plotFlag && mod(i,10)==0
            vv = Vp*p;
            clf; color_line3(xp,yp,vv,vv,'.');
            axis equal; axis tight; colorbar
            title(sprintf('time = %f, max solution val = %f',i*dt,max(abs(vv(:)))))
            drawnow
        end
        if (mod(i,ceil(Nstep/5))==0)
            disp(sprintf('backwards prop: tstep %d out of %d',i,Nstep))
        end        
    end
    
    
    if useTraces
        p = p_old - p; % w = (I-A*L)*p0 = p0 - A*L*p0
        p_recon = p_recon + p; % accumulate terms!
    else        
        p = p0 + p;
        p_recon = p;
    end
    
    % plot initial cond
    p0ex = Vq*Pq*pex(xq,yq); % use projection instead of exact?
    err = diag(wq)*(Vq*J).*(Vq*(p_recon)-p0ex).^2;
    L2err(iter+1) = sqrt(sum(err(:)));
    uq = diag(wq)*(Vq*J).*pex(xq,yq).^2;
    Unorm = sqrt(sum(uq(:)));
    
    disp(['iter ' num2str(iter) ', L2 err in recon = ',num2str(L2err/Unorm)])
    %     pause
end

keyboard

%%

figure
vv = pex(xp,yp)-Vp*p_recon;
% vv = Vp*p_recon;
clf; color_line3(xp,yp,vv,vv,'.');
axis equal; axis tight; colorbar
title(sprintf('iter = %d, reconstruction error = %f',iter,L2err(iter)))

%%

return


function [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,dt,bcOpt,pbc)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

global tau0
tau = sign(dt)*tau0;

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

% evaluate upwind fluxes
ndotdU = nx.*du + ny.*dv;

if strcmp(bcOpt,'d') % Impose reflective boundary conditions (p+ = -p-)
    ndotdU(mapB) = 0;
    dp(mapB) = 2*(-p(vmapB));
elseif strcmp(bcOpt,'n') % free surface BCs
    ndotdU(mapB) = -2*(nx(mapB).*u(vmapM(mapB)) + ny(mapB).*v(vmapM(mapB)));
    dp(mapB) = 0;
elseif strcmp(bcOpt,'abc') % basic abcs
    un = nx(mapB).*u(vmapM(mapB)) + ny(mapB).*v(vmapM(mapB));
    dp(mapB) = -p(vmapB);
    ndotdU(mapB) = -un;
end

if strcmp(bcOpt,'rd')
    dp(mapB) = 2*(pbc-p(vmapB));
    ndotdU(mapB) = 0;
end

fluxp =  tau*dp - ndotdU;
fluxu =  (tau*ndotdU - dp);


pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
rhsp =  -divU + LIFT*(Fscale.*fluxp)/2.0;
rhsu =  -dpdx + LIFT*(Fscale.*fluxu.*nx)/2.0;
rhsv =  -dpdy + LIFT*(Fscale.*fluxu.*ny)/2.0;

global Pq cq Vq
rhsp = Pq*(cq.*(Vq*rhsp));


return;

