clear

% Order of polymomials used for approximation
N = 5;
FinalTime = 2;
r = 1;
CFL = 2^(-r);

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

%% make 2-element mesh 

xhat = [-1 + .5*(1+r) .5*(1+r)];
K = size(xhat,2);
Jhat = .5;

%% moving mapping

a = 0.25;
Phi = @(x,t) x + a*sin(pi*t).*(1-x).*(1+x);
dPhi_t = @(x,t) a*pi*cos(pi*t).*(1-x).*(1+x);

% xhat = r;
% K = 1;

%% time stepping

% Set initial conditions
u = 1+exp(-10*xhat.^2);
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
xmin = min(min(diff(xhat,1)));
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
        J = Dr*Phi(xhat,time_local); 
        x = Phi(xhat,time_local);        
        xt = dPhi_t(xhat,time_local);
        Jrt = xt./Jhat; % chain rule: J*dxhat/dt = J*dxhat/dx * dx/dt = dx/dt since J*dxhat/dx = 1 in 1D
                
        Jruf = (Vf*Jrt).*(Vf*u);
        mapP = [4 3; 2 1]';
        flux = Jruf(mapP)-Jruf;
        Jxth = Vh*Jrt;
        uh = Vh*u;
        rhsuJ = Vh'*(.5*Qh*(Jxth.*uh) + .5*Jxth.*(Qh*uh) + .5*(Qh*Jxth).*uh) + .5*Vf'*([-1;1].*flux);
        rhsuJ = M\rhsuJ;        
        
        rhsJ = Dr*Jrt;
                
        rhstest = sum(sum(u.*(M*rhsuJ) - (Pq*(.5*(Vq*u).^2)).*(M*rhsJ)));
        
        resJ  = rk4a(INTRK)*resJ + dt*rhsJ;
        resuJ = rk4a(INTRK)*resuJ + dt*rhsuJ;
        
        uJ = uJ + rk4b(INTRK)*resuJ;
        Jh = Jh + rk4b(INTRK)*resJ;
        
        for e = 1:K            
            MJ = Vq'*diag(wq.*(Vq*Jh(:,e)))*Vq;
            u(:,e) = MJ\(M*uJ(:,e));
        end

%         u = Pq*((Vq*uJ)./(Vq*Jh)); 
        
    end  
        
    unorm(tstep) = sum(sum(diag(wq)*((Vq*Jh).*(Vq*u).^2)));
    
    val = 0;
    for e = 1:K
        MinvJ = Vq'*diag(wq./(Vq*Jh(:,e)))*Vq;
        val = val + uJ(:,e)'*MinvJ*uJ(:,e);
    end
    unorm2(tstep) = val;
    
    if  mod(tstep,5)==0 || tstep==Nsteps
        time = (tstep)*dt;
        plot(Vp*x,Vp*u,'-');
%         plot(Vp*r,Vp*Jh,'--');
%         plot(Vp*x,Vp*u,'-');
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



