function [u] = Advec3D(u, FinalTime)

% function [u] = Advec3D(u, FinalTime)
% Purpose  : Integrate 3D advection until FinalTime starting with initial condition, u

Globals3D;
time = 0;

% Runge-Kutta residual storage  
resu = zeros(Np,K); 

% compute time step size
dt = .125/max( ((N+1)^2)*.5*Fscale(:));

% outer time step loop 
tstep = 1;
while (time<FinalTime)

  if(time+dt>FinalTime), dt = FinalTime-time;  end;

  % low storage RK
  for INTRK = 1:5
    timelocal = time + rk4c(INTRK)*dt;
    rhsu = AdvecRHS3D(u,timelocal);
    resu = rk4a(INTRK)*resu + dt*rhsu;
    u = u+rk4b(INTRK)*resu;
  end;

  % Increment time
  time = time+dt; tstep = tstep+1;
end;
return

function [rhsu] = AdvecRHS3D(u,time)

% function [rhsu] = AdvecRHS3D(u,time)
% Purpose  : Evaluate RHS flux in 3D advection

Globals3D;

% form field differences at faces
alpha=1; du = zeros(Nfp*Nfaces,K); 
du(:) = 0.5*(u(vmapM)-u(vmapP)).*(nx(:) - alpha*abs(nx(:)));

% impose boundary condition at x=0
ubc  = 0;%exp(-1*( (Fx(mapB)-time).^2 + Fy(mapB).^2 + Fz(mapB).^2));
du(mapB) = 0.5*(u(vmapB)-ubc).*(nx(mapB)-alpha*abs(nx(mapB)));

% compute right hand sides of the semi-discrete PDE
rhsu = -(rx.*(Dr*u)+sx.*(Ds*u)+tx.*(Dt*u)) + LIFT*(Fscale.*(du));

return