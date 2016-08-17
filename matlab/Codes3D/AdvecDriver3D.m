% Driver script for solving the 3D advection equations
Globals3D;

% Order of polymomials used for approximation 
N = 1;

% Generate simple mesh
filename = 'Grid/cube1.neu';
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D(filename);

% Initialize solver and construct grid and metric
StartUp3D;

% for plotting - build coordinates of all the nodes
[rp sp tp] = EquiNodes3D(20);
va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)'; vd = EToV(:,4)';
xp = 0.5*(-(1+rp+sp+tp)*VX(va)+(1+rp)*VX(vb)+(1+sp)*VX(vc)+(1+tp)*VX(vd));
yp = 0.5*(-(1+rp+sp+tp)*VY(va)+(1+rp)*VY(vb)+(1+sp)*VY(vc)+(1+tp)*VY(vd));
zp = 0.5*(-(1+rp+sp+tp)*VZ(va)+(1+rp)*VZ(vb)+(1+sp)*VZ(vc)+(1+tp)*VZ(vd));
ids = yp > 0; 
xp = xp(ids); yp = yp(ids); zp = zp(ids);
Ip = Vandermonde3D(N,rp,sp,tp)*invV;

% set initial conditions
u = exp(-2.5*(x.^2 + y.^2 + z.^2));

% Solve Problem
FinalTime = 3;
% Runge-Kutta residual storage  
resu = zeros(Np,K); 

% compute time step size
dt = .125/max( ((N+1)^2)*.5*Fscale(:));

% outer time step loop 
time = 0;
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
  
  if (mod(tstep,50)==0)
      up = Ip*reshape(u,Np,K);
      color_line3(xp,yp,zp,up(ids),'.');
      title(sprintf('time t = %f',time))
      view(3); 
      caxis([0,1])
      drawnow
  end

  % Increment time
  time = time+dt; tstep = tstep+1;
end
