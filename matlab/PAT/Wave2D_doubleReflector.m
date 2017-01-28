function Wave2D_simple

Globals2D

% generate triangle
yt = 13.243323833946539;
VX = [14.9071 35.3 -10]; 
VY = [0 yt yt];
K = 1; EToV = 1:3;
N = 5;

StartUp2D

BCType = ones(size(EToV));
nref = 3;
for ref = 1:nref
    Refine2D(ones(size(EToV)));
    StartUp2D
end
% PlotMesh2D;return

% fine-mesh nodes for plotting
[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

%% pick out top faces

global vmapB_top mapB_top

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

% plot(x,y,'.')
% hold on
% xf = x(Fmask(:),:); yf = y(Fmask(:),:);
% plot(x(vmapB_top),y(vmapB_top),'ro')
% plot(xf(mapB_top),yf(mapB_top),'ro')

%% set up initial condition

x0 = 14; y0 = 7;
p = exp(-1^2*((x-x0).^2 + (y-y0).^2));
u = zeros(Np, K); 
v = zeros(Np, K);

p = exp(-2*(y-y0).^2).*(abs(x-x0)<8);

%% set up time integration parameters

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

global c2 
c2 = 1.5^2; % wave speed c = 1.5?

% compute time step size
CN = (N+1)*(N+2)/2; % trace constant on the reference element
CNh = max(c2*CN*max(Fscale(:)));
dt = 2/CNh; % cfl condition of 2

%% run simulation 

FinalTime = 12.20;

Nstep = ceil(FinalTime/dt)
dt = FinalTime/Nstep

trace = zeros(length(vmapB_top),Nstep);

% outer time step loop
for tstep = 1:Nstep

    time = tstep*dt;
    for INTRK = 1:5
        
        timelocal = time + rk4c(INTRK)*dt;
        
        [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,timelocal);
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        p = p+rk4b(INTRK)*resp;
                
    end;
    
    if mod(tstep,10)==0
        clf
        vv = Vp*p;
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f',time))
        drawnow
    end        
    
    trace(:,tstep) = p(vmapB_top);
end

% %%
% figure
% hold on
% for i = 1:Nstep;    
%     vv = trace(:,i);
%     color_line3(x(vmapB_top),i*dt*ones(size(vv)),vv,vv,'.');
% end
%%
% keyboard



function [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

% evaluate upwind fluxes
ndotdU = nx.*du + ny.*dv;

% Impose reflective boundary conditions (neumann)
dp(mapB) = 0;
ndotdU(mapB) = -2*(nx(mapB).*u(vmapB)+ny(mapB).*v(vmapB));

global vmapB_top mapB_top
dp(mapB_top) = -p(vmapB_top);
ndotdU(mapB_top) = -(nx(mapB_top).*u(vmapB_top)+ny(mapB_top).*v(vmapB_top));

% compute numerical fluxes
tau = 1;
fluxp =  tau*dp - ndotdU;
fluxu =  (tau*ndotdU - dp).*nx;
fluxv =  (tau*ndotdU - dp).*ny;

% compute derivatives
pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

rhsp =  -divU + LIFT*(Fscale.*fluxp)/2.0;
rhsu =  -dpdx + LIFT*(Fscale.*fluxu)/2.0;
rhsv =  -dpdy + LIFT*(Fscale.*fluxv)/2.0;

global c2
rhsp = c2*rhsp;

return;

