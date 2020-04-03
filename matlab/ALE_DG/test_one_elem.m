clear

% Order of polymomials used for approximation
N = 15;
FinalTime = 2;
CFL = .125;

r = JacobiGL(0,0,N);
V = Vandermonde1D(N,r);
Dr = GradVandermonde1D(N,r)/V;

[rq wq] = JacobiGQ(0,0,N+1);
% [rq wq] = JacobiGL(0,0,N);
Vq = Vandermonde1D(N,rq)/V;
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq));

rf = [-1;1];
Vf = Vandermonde1D(N,rf)/V;

rp = linspace(-1,1,500)';
Vp = Vandermonde1D(N,rp)/V;

Qr = Pq'*M*Dr*Pq;
E = Vf*Pq;
B = diag([-1;1]);

Qh = .5*[Qr-Qr' E'*B;
    -B*E B];
Vh = [Vq;Vf];

%% moving mapping

a = 0.25;
Phi = @(x,t) x + a*sin(pi*t).*(1-x).*(1+x);
dPhi_t = @(x,t) a*pi*cos(pi*t).*(1-x).*(1+x);

x = r;

%% time stepping

% Set initial conditions
u = 1+exp(-10*r.^2);
% u = exp(sin(pi*x));
uJ = u;
Jh = ones(size(u));

rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0];

% Runge-Kutta residual storage
resuJ = zeros(size(u));
resJ = zeros(size(u));

% compute time step size
xmin = min(diff(r,1));
dt  = CFL*xmin;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

time_local = 0;

figure(1)
for tstep=1:Nsteps
    
    time = (tstep-1)*dt;
    
    for INTRK = 1:5
        
        time_local = time + rk4c(INTRK)*dt;
        
        % d(uJ)/dt + div(J*xt*u)
        J = Dr*Phi(r,time_local); 
        x = Phi(r,time_local);        
        xt = dPhi_t(r,time_local);
        Jrt = xt; % chain rule: J*dxhat/dt = J*dxhat/dx * dx/dt = dx/dt since J*dxhat/dx = 1 in 1D
                
        flux = (Vf*Jrt).*(Vf*u);
        Jxth = Vh*Jrt;
        uh = Vh*u;
%         uh2 = Vh*Pq*((Vq*uJ)./(Vq*J));
        uh2 = uh;
        rhsuJ = Vh'*(.5*Qh*(Jxth.*uh) + .5*Jxth.*(Qh*uh) + .5*(Qh*Jxth).*uh2) + .5*Vf'*([-1;1].*flux);
        rhsuJ = M\rhsuJ;
        
        rhsJ = Dr*Jrt;
                
        rhstest = sum(u.*(M*rhsuJ) - (Pq*(.5*(Vq*u).^2)).*(M*rhsJ));
        
        resJ  = rk4a(INTRK)*resJ + dt*rhsJ;
        resuJ = rk4a(INTRK)*resuJ + dt*rhsuJ;
        
        uJ = uJ + rk4b(INTRK)*resuJ;
        Jh = Jh + rk4b(INTRK)*resJ;
        
        MJ = Vq'*diag(wq.*(Vq*Jh))*Vq;
        u = MJ\(M*uJ);

%         u = Pq*((Vq*uJ)./(Vq*Jh)); 
        
    end  
    
    MJ = Vq'*diag(wq.*(Vq*Jh))*Vq;
    unorm(tstep) = u'*MJ*u;
    MinvJ = Vq'*diag(wq./(Vq*Jh))*Vq;
    unorm2(tstep) = uJ'*MinvJ*uJ;
    
    if  mod(tstep,25)==0 || tstep==Nsteps
        time = (tstep)*dt;
        plot(Vp*x,Vp*u,'-');
        hold on
        plot(x,u,'o');
        axis([-1,1,-1.1,2.1])
        hold off        
        title(sprintf('time = %g, rhstest = %g\n',tstep*dt,rhstest))
        drawnow
    end
    
end

figure(2)
semilogy(dt*(1:Nsteps),abs(unorm-unorm(1)),'--')
hold on
semilogy(dt*(1:Nsteps),abs(unorm2-unorm2(1)),'-.')



