clear
Globals2D

kd = [4 8 16 32]% 64]
h= 2./kd

N = 4;
M = 1;
k = 1;
for i = 1:length(kd)
    
    K1D = kd(i);
    FinalTime = .3;
    [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
    
    iids = abs(VX)<1 & abs(VY) < 1;
    VX(iids) = VX(iids) + .125/K1D*randn(size(VX(iids)));
    VY(iids) = VY(iids) + .125/K1D*randn(size(VX(iids)));
    
    StartUp2D;
%     PlotMesh2D;return
    
    
    % Set up wavespeed function
    %cfun = @(x,y) ones(size(x));
    cfun = @(x,y) 1+ 0.5*sin(pi*x).*sin(pi*y); % smooth velocity
    %cfun = @(x,y) (1 + .5*sin(2*pi*x).*sin(2*pi*y) + (y > 0)); % piecewise smooth velocity
    
    % Set periodic manufactured solution
    pfun = @(x,y,t) sin(k*pi*x).*sin(k*pi*y).*cos(k*pi*t);
    ufun = @(x,y,t) -cos(k*pi*x).*sin(k*pi*y).*sin(k*pi*t);
    vfun = @(x,y,t) -sin(k*pi*x).*cos(k*pi*y).*sin(k*pi*t);
    ffun = @(x,y,t)  k*pi*(-1./(cfun(x,y))+2).*sin(k*pi*x).*sin(k*pi*y).*sin(k*pi*t);
    
    % generate quadrature points
    global Pq CqM Vq
    Nq = 2*N+1;
    [rq sq wq] = Cubature2D(Nq); % integrate u*v*c
    Vq = Vandermonde2D(N,rq,sq)/V;
    xq = Vq*x; yq = Vq*y;
    
    % construct the projection matrix for nodal basis Pq
    Pq = V*V'*Vq'*diag(wq);
    
    % construct matrix Cq
    Cq = cfun(xq,yq);
    
    % construct the projection matrix for cfun to degree M
    VMq = Vandermonde2D(M,rq,sq);
    CqM = VMq*VMq'*diag(wq)*Cq;
%     CqM = Cq;
    
    % initial condition
    p = pfun(x,y,0);
    u = ufun(x,y,0);
    v = vfun(x,y,0);
    
    time = 0;
    
    % Runge-Kutta residual storage
    resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);
    
    % compute time step size
    CN = (N+1)*(N+2)/2; % trace inequality constant
    CNh = max(CN*max(Fscale(:)));
    dt = 2/CNh;
    
    % outer time step loop
    tstep = 0;
    
    while (time<FinalTime)
        if(time+dt>FinalTime), dt = FinalTime-time; end
        
        for INTRK = 1:5
            
            timelocal = time + rk4c(INTRK)*dt;
            
            f = Pq*ffun(xq,yq,timelocal);            
            
            [rhsp, rhsu, rhsv] = acousticsRHS2D_compare(p,u,v,f);
            
            % initiate and increment Runge-Kutta residuals
            resp = rk4a(INTRK)*resp + dt*rhsp;
            resu = rk4a(INTRK)*resu + dt*rhsu;
            resv = rk4a(INTRK)*resv + dt*rhsv;
            
            % update fields
            u = u+rk4b(INTRK)*resu;
            v = v+rk4b(INTRK)*resv;
            p = p+rk4b(INTRK)*resp;
            
        end
        
        % Increment time
        time = time+dt; tstep = tstep+1;
        if mod(tstep,10)==0            
            fprintf('on timestep %d out of %d\n',tstep,ceil(FinalTime/dt))
        end
    end
    
    p_exact = pfun(x,y,FinalTime);
    p_exact_quadrature = pfun(xq,yq,FinalTime);
    
    p_quadrature = Vq * p;
    
    [d1,d2]=size(p_quadrature);
    error_accumulation = 0;
    for j1=1:d2
        for j2=1:d1
            err = p_quadrature(j2,j1)-p_exact_quadrature(j2,j1);
            error_accumulation = error_accumulation + err*err*wq(j2)*J(1,j1);
        end
    end
    
    error_l2(i) = sqrt(error_accumulation);
    error_fro(i) = norm(p-p_exact,'fro');
end



function [rhsp, rhsu, rhsv] = acousticsRHS2D_compare(p,u,v,f)

Globals2D;

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

% evaluate upwind fluxes
ndotdU = nx.*du + ny.*dv;

% Impose reflective boundary conditions (p+ = -p-)
ndotdU(mapB) = 0;
dp(mapB) = -2*p(vmapB);

tau = 1;
fluxp =  tau*dp - ndotdU;
fluxu =  (tau*ndotdU - dp).*nx;
fluxv =  (tau*ndotdU - dp).*ny;

pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
rhsp =  -divU + f + LIFT*(Fscale.*fluxp)/2.0;
rhsu =  -dpdx + LIFT*(Fscale.*fluxu)/2.0;
rhsv =  -dpdy + LIFT*(Fscale.*fluxv)/2.0;

global Pq CqM Vq
rhsp = Pq*(CqM.*(Vq*rhsp));
return;
end