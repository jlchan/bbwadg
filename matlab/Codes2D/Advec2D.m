function Advec2D

Globals2D

% Polynomial order used for approximation 
N = 3;

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/squarereg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/squareireg.neu');
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/block2.neu');

% Initialize solver and construct grid and metric
StartUp2D;

% rebuild maps for periodic
BuildPeriodicMaps2D(2,2);

% % Set initial conditions
cx = .1; cy = .1;
D = (x-cx).^2 + (y-cy).^2;
u = exp(-D*5^2).*(1-x.^2).*(1-y.^2);

% run 
FinalTime = 2;
[u,time] = AdvecRun2D(u,FinalTime);


function [u,time] = AdvecRun2D(u, FinalTime)

Globals2D;
time = 0;

% Runge-Kutta residual storage  
resu = zeros(Np,K); 

% compute time step size
rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D; dt = min(dtscale)*rmin*2/3

% outer time step loop 
while (time<FinalTime)
  
  if(time+dt>FinalTime), dt = FinalTime-time; end

   for INTRK = 1:5    
      % compute right hand side of TM-mode Maxwell's equations
      alpha = 1; du = zeros(Nfp*Nfaces,K);
      du(:) = 0.5*(u(vmapM)-u(vmapP)).*(nx(:) - alpha*abs(nx(:)));
      rhsu = -(rx.*(Dr*u)+sx.*(Ds*u)) + LIFT*(Fscale.*(du));

      % initiate and increment Runge-Kutta residuals
      resu = rk4a(INTRK)*resu + dt*rhsu;  
        
      % update fields
      u = u+rk4b(INTRK)*resu; 
   end
   
   clf;color_line3(x,y,u,u,'.');drawnow
   
   % Increment time
   time = time+dt;
end
return

