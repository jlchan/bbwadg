Globals2D

N = 4;
K1D = 8;
FinalTime = 10;

[Nv, VX, VY, K, EToV] = QuadMesh2D(K1D);
VX = (1+VX)/2; VY = (1+VY)/2;

StartUp2D;

BuildPeriodicMapsQuad2D(1,1);

% plotting nodes
[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V; 
xp = Vp*x; yp = Vp*y;

LIFT = LiftQuad2D_DGSEM();

%%

% advection equation: du/dt + bx * du/dx + by * du/dy = 0
global bx by

opt = 2; 
if opt==1 % warped flow on [-1,1]^2
    
    % beta = (bx,by)
    bx = 1+0*cos(pi*y);
    by = 0*sin(pi*x);
    u = exp(-50*((x-.5).^2+(y-.5).^2)); % Gaussian pulse
    
elseif opt==2 % swirling flow ("double glazing") on [0,1]^2
    
    bx = pi*sin(pi*x).*cos(pi*y);
    by = -pi*cos(pi*x).*sin(pi*y);
    u = exp(x.*y);
    
end

%%


time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
CNh = max(CN*max(Fscale(:)));
dt = 1/CNh;

% outer time step loop
figure
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;
for i = 1:Nsteps    
    for INTRK = 1:5
        
        timelocal = time + rk4c(INTRK)*dt;
        
        rhsu  = advecRHS2D(u,timelocal);
        
        % initiate and increment Runge-Kutta residuals
        resu = rk4a(INTRK)*resu + dt*rhsu;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        
    end;
    
    if mod(i,10)==0 || i==Nsteps
        clf
        pp = u;
        vv = Vp*pp;        
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f',dt*i))
        drawnow
    end
        
end


function [rhsu] = advecRHS2D(u,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;
global bx by

tau = 1;
uP = u(vmapP); uM = u(vmapM); 
ujump = uP-uM;

bxf = bx(vmapM);
byf = by(vmapM);
bn = bxf.*nx(:) + byf.*ny(:);

fx = bx.*u; fy = by.*u;
dfx = rx.*(Dr*fx) + sx.*(Ds*fx);
dfy = ry.*(Dr*fy) + sy.*(Ds*fy);
divF = dfx + dfy;

ur = (Dr*u); us = (Ds*u);
dudx = rx.*ur + sx.*us;
dudy = ry.*us + sy.*us;
bGradu = bx.*dudx + by.*dudy;

% simplifications: bx, by = continuous, combine flux terms
fluxu = .5*bn.*ujump - tau*.5*abs(bn).*ujump;
fluxu = reshape(fluxu,Nfp*Nfaces,K);
rhsu =  -(.5*(divF+bGradu) + LIFT*(Fscale.*fluxu));

end

