function Wave2D

% clear all, clear
clear -global *

Globals2D

N = 4;
global Pq cq Vq c2
c2 = 1.5^2;
cfun = @(x,y) c2*ones(size(x));
% cfun = @(x,y) 1 + (x > 0);
% cfun = @(x,y) 1 + .25*sin(2*pi*x).*sin(2*pi*y); % smooth velocity
% cfun = @(x,y) (1 + .25*sin(2*pi*x).*sin(2*pi*y) + (y > 0)); % piecewise smooth velocity

FinalTime = 12.2;
% FinalTime = 15;

% generate triangle
nref = 6;
yt = 13.243323833946539;
% yt = 12.5 - (yt - 12.5)
VX = [14.9071 35.3 -10]; VY = [0 yt yt]; 
% VX = [14.9071,55.6929,-34.9070]; VY = [-yt,yt,yt]; % bigger tri?
K = 1; EToV = 1:3;
StartUp2D
BCType = ones(size(EToV));
for ref = 1:nref
    Refine2D(ones(size(EToV)));
    StartUp2D
end

% PlotMesh2D; axis on;return

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

cq = cfun(xq,yq);

% FinalTime = 2/min(cq(:));

%% pick out top faces

global vmapB_top mapB_top vmapB_bot mapB_bot

vmapB_top = [];
mapB_top = [];
for e = 1:K
    for f = 1:Nfaces
        yf = y(Fmask(:,f),e);
        if norm(yf-yt) < 1e-6
            vids = Fmask(:,f) + Np*(e-1); % volume id
            fids = (1:Nfp) + Nfp*(f-1) + Nfp*Nfaces*(e-1); % face id
            mapB_top = [mapB_top; fids(:)];
            vmapB_top = [vmapB_top; vids(:)];
        end
    end
end

mapB_bot = setdiff(mapB,mapB_top);
vmapB_bot = vmapM(mapB_bot);

% make top boundary quadrature points
V1D = Vandermonde1D(N,JacobiGL(0,0,N));
[r1Dq w1Dq] = JacobiGQ(0,0,N);
Vq1D = Vandermonde1D(N,r1Dq)/V1D;

Pq1D = V1D*V1D'*Vq1D'*diag(w1Dq);
xq_top = Vq1D*x(reshape(vmapB_top,Nfp,length(vmapB_top(:))/Nfp));
xtop = Pq1D*xq_top;
% plot(x(vmapB_bot),y(vmapB_bot),'o');return

%% params setup

x0 = 14; y0 = 7;

% gaussian pulse
pex = @(x,y) exp(-1^2*((x-x0).^2 + (y-y0).^2));
% pex = @(x,y) cos(.5*pi*x).*cos(.5*pi*y);

% smooth annulus
sig = @(x) 1-1./(1+exp(-1e2*x));
r2 = @(x,y) x.^2 + y.^2;
pex = @(x,y) sig((.025*r2(x-x0,y-y0)-.25).^2);
% pex = @(x,y) abs(y-5.5)<.2; 

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
dt = 2/CNh;

Nstep = ceil(FinalTime/dt);
dt = FinalTime/Nstep;

tau = 1;

plotFlag = 1;

%% generate synthetic data

NstepRK = Nstep*5;
Ubc{1} = zeros(length(mapB_top(:)),NstepRK);
Ubc{2} = zeros(length(mapB_top(:)),NstepRK);

if 1
    for i = 1:Nstep
        for INTRK = 1:5
            [rhsp, rhsu, rhsv,UbcI] = acousticsRHS2D(p,u,v,tau,'abc');
            resp = rk4a(INTRK)*resp + dt*rhsp;
            resu = rk4a(INTRK)*resu + dt*rhsu;
            resv = rk4a(INTRK)*resv + dt*rhsv;
            p = p+rk4b(INTRK)*resp;
            u = u+rk4b(INTRK)*resu;
            v = v+rk4b(INTRK)*resv;
            
            % save traces
            id0 = ((i-1)*5+(INTRK-1)); % zero index
            a = 1e-1; % noise threshhold
            Ubc{1}(:,NstepRK - id0) = UbcI{1} + a*randn(size(UbcI{1}));
            Ubc{2}(:,NstepRK - id0) = UbcI{2} + a*randn(size(UbcI{1}));
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
    
else
    
    tRK = (0:Nstep-1)*dt;
    tRK = repmat(dt*rk4c(:),1,Nstep) + repmat(tRK,5,1);
    tRK = tRK(:);
    load PAT/TripleHair_Param.mat
    load PAT/TripleHair_Data.mat

    dt_data = 12.2 / size(PA_data,1);
    xn = xn + x0_shift;
    tt = 0:dt_data:12.2-dt_data;
    
    % extend data by 0 outside interval for now?
    xx = linspace(-10,35,(35-(-10))/.125+1);
    %     data = zeros(size(PA_data,1),length(xx));
    [~,startid] = min(abs(xx-min(xn(:))));
    [~,endid] = min(abs(xx-max(xn(:))));
    ids = startid:endid;
    
    if 1 % try interpolating in time first        
        data_t = zeros(NstepRK,length(xx));
        for i = 1:size(PA_data,2)
            tvec_pt = csapi(tt,PA_data(:,i),tRK);
            
%             fc = 20; fs = 80;
%             [b,a] = butter(6,fc/(fs/2));
%             tvec_pt = filter(b,a,tvec_pt);
            data_t(:,ids(i)) = tvec_pt;
            if mod(i,round(size(PA_data,1)/10))==0
                disp(sprintf('interp in time: on i = %d out of %d\n',i,size(data_t,2)));
            end
        end        
        
        for i = 1:NstepRK                        
            yq_top = csapi(xx,data_t(NstepRK-i+1,:),xq_top); % flip to read from end (last time) first
            tmp = Pq1D*yq_top;
            Ubc{1}(:,i) = tmp(:);
            if mod(i,round(size(PA_data,1)/10))==0
                disp(sprintf('interp in space: on i = %d out of %d\n',i,NstepRK));
            end
        end
    else
        disp('Interpolating data in space')
        data1 = zeros(size(PA_data,1),length(mapB_top(:)));
        
        % interpolate data in space first
        yy = zeros(size(xx));
        for i = 1:size(PA_data,1)
            yy(ids) = PA_data(i,:); % flip data
            yq_top = csapi(xx,yy,xq_top);
            tmp = Pq1D*yq_top;
            data1(i,:) = tmp(:);
            if mod(i,round(size(PA_data,1)/10))==0
                disp(sprintf('interp in space: on i = %d out of %d\n',i,size(PA_data,1)));
            end
        end
        
        disp('Interpolating data in time')
        for i = 1:size(data1,2)
            tvec_pt = csapi(tt,data1(:,i),tRK);
            Ubc{1}(i,:) = tvec_pt;
            if mod(i,round(size(PA_data,1)/10))==0
                disp(sprintf('interp in time: on i = %d out of %d\n',i,size(data1,2)));
            end
        end
    end
%     for i = 1:size(Ubc{1},2)
%         Ubc{1}(1:startid-1,i) = Ubc{1}(startid,i)*ones(startid-1,1);
%         keyboard
%         Ubc{1}(endid+1:size(Ubc{1},2),i) = Ubc{1}(endid,i)*ones(size(Ubc{1},2)-endid-1,1);
%     end
%         for i = 1:NstepRK
%             clf
%             plot(x(vmapB_top),Ubc{1}(:,i),'o')
% %                     clf
% %             plot(x(vmapB_top),tmp(:),'o')
% %             hold on
% %             plot(xx,yy,'-')
%             axis([-10 35 -550 550])
%             pause(.01)
%         end
%     Ubc{1} = Ubc{1}/ max(abs(Ubc{1}(:))); % normalize
    
    keyboard
    
end


%% playback : reverse propagate to generate Rm from pressure boundary data

dt = -dt;
p = zeros(Np,K); u = p; v = p;
resp = zeros(Np,K); resu = resp; resv = resp;

for i = 1:Nstep
    for INTRK = 1:5
        id0 = ((i-1)*5+(INTRK-1)); % 1-index
        UbcI{1} = Ubc{1}(:,id0+1);
        UbcI{2} = Ubc{2}(:,id0+1);
        
        [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,sign(dt)*tau,'rd',UbcI);
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        p = p + rk4b(INTRK)*resp;
        u = u + rk4b(INTRK)*resu;
        v = v + rk4b(INTRK)*resv;
    end;
    
    if (mod(i,ceil(Nstep/5))==0)
        disp(sprintf('first backwards prop: tstep %d out of %d',i,Nstep))
    end
    if plotFlag && mod(i,10)==0
        vv = Vp*p;
        clf; color_line3(xp,yp,vv,vv,'.');
        axis equal; axis tight; colorbar
        caxis([-2.5 2.5])
        title(sprintf('first reversal time = %f, max solution val = %f',i*dt,max(abs(vv(:)))))
        drawnow
    end
end

vv = Vp*p;
clf; color_line3(xp,yp,vv,vv,'.');
axis equal; axis tight; colorbar
title(sprintf('first reversal time = %f, max solution val = %f',i*dt,max(abs(vv(:)))))
keyboard

% save reversed data
Rm = p;


%% iterate

% initalize pressure condition
p0 = Rm;

for iter = 1:10
    
    p = p0; u = zeros(Np,K); v = zeros(Np,K);
    resp = zeros(Np,K); resu = resp; resv = resp;
    
    % compute forward propagation
    dt = abs(dt);
    for i = 1:Nstep
        for INTRK = 1:5
            [rhsp, rhsu, rhsv, UbcI] = acousticsRHS2D(p,u,v,sign(dt)*tau,'abc');
            resp = rk4a(INTRK)*resp + dt*rhsp;
            resu = rk4a(INTRK)*resu + dt*rhsu;
            resv = rk4a(INTRK)*resv + dt*rhsv;
            p = p + rk4b(INTRK)*resp;
            u = u + rk4b(INTRK)*resu;
            v = v + rk4b(INTRK)*resv;
        end;
        
        if (mod(i,ceil(Nstep/5))==0)
            disp(sprintf('fwd prop: tstep %d out of %d',i,Nstep))
        end
        %         if mod(i,10)==0
        %             vv = Vp*p;
        %             clf; color_line3(xp,yp,vv,vv,'.');
        %             title(sprintf('forwards time = %f, max solution val = %f',time,max(abs(vv(:)))))
        %             axis equal; axis tight; colorbar
        %             drawnow
        %         end
    end
    
    %     % compute new auxiliary initial conditions
    %     p = 0*p;
    %     u = 0*u;
    %     v = 0*v;
    
    % compute backwards propagation
    dt = -dt;
    for i = 1:Nstep
        for INTRK = 1:5
            [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,sign(dt)*tau,'d');
            
            resp = rk4a(INTRK)*resp + dt*rhsp;
            resu = rk4a(INTRK)*resu + dt*rhsu;
            resv = rk4a(INTRK)*resv + dt*rhsv;
            p = p + rk4b(INTRK)*resp;
            u = u + rk4b(INTRK)*resu;
            v = v + rk4b(INTRK)*resv;
        end;
        
        if (mod(i,ceil(Nstep/5))==0)
            disp(sprintf('backwards prop: tstep %d out of %d',i,Nstep))
        end
        %         if mod(i,10)==0
        %             vv = Vp*p;
        %             clf; color_line3(xp,yp,vv,vv,'.');
        %             title(sprintf('backwards, time = %f, max solution val = %f',time,max(abs(vv(:)))))
        %             axis equal; axis tight; colorbar
        %             drawnow
        %         end
    end
    
    % update initial condition
    %     p0 = (p0-p) + Rm;
    p0 = p + Rm;
    
    % plot initial cond
    err = diag(wq)*(Vq*J).*(Vq*p0-pex(xq,yq)).^2;
    L2err(iter) = sqrt(sum(err(:)));
    uq = diag(wq)*(Vq*J).*pex(xq,yq).^2;
    Unorm = sqrt(sum(uq(:)));
    
    figure
    vv = Vp*p0;
    clf; color_line3(xp,yp,vv,vv,'.');
    axis equal; axis tight; colorbar
    title(sprintf('iter = %d, reconstruction error = %f',iter,L2err(iter)))
    
    disp(['iter ' num2str(iter) ', L2 err in recon = ',num2str(L2err/Unorm)])
%     keyboard
end

keyboard
return


function [rhsp, rhsu, rhsv, UbcO] = acousticsRHS2D(p,u,v,tau,bcOpt,Ubc)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

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

global vmapB_top mapB_top vmapB_bot mapB_bot

% impose pressure BC
if strcmp(bcOpt,'rd')
    dp(mapB_top) = 2*(Ubc{1}-p(vmapB_top));
    ndotdU(mapB_top) = 0;
end
if nargout > 3
    UbcO{1} = p(vmapB_top);
    UbcO{2} = ndotdU(mapB_top);
end

% % dirichlet
% dp(mapB_bot) = -2*p(vmapB_bot);
% ndotdU(mapB_bot) = 0;

% neumann for reflectors
dp(mapB_bot) = 0;
ndotdU(mapB_bot) = -2*(nx(mapB_bot).*u(vmapB_bot) + ny(mapB_bot).*v(vmapB_bot));

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

global Pq cq Vq c2
% rhsp = Pq*(cq.*(Vq*rhsp));
rhsp = c2*rhsp;

return;

