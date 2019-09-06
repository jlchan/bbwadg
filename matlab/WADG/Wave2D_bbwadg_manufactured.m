clear
Globals2D

kd = [4 8 16 32];
h = 2./kd;
N = 1;
MM = 4;
CFL = .750;

k = 1;
for i = 1:length(kd)
    
    K1D = kd(i);
    FinalTime = .7;
    [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
    StartUp2D;
%     BuildPeriodicMaps2D(2,2);
    
    %% Set up wavespeed function
    %cfun = @(x,y) ones(size(x));
    cfun = @(x,y) 1 + 0.5*sin(1+pi*x).*sin(2+.5*pi*y); % smooth velocity
    %     cfun = @(x,y) (1 + .5*sin(2*pi*x).*sin(2*pi*y) + (y > 0)); % piecewise smooth velocity
    
    %% Set periodic manufactured solution
    pfun = @(x,y,t)  sin(k*pi*x).*sin(k*pi*y).*cos(k*pi*t);
    ufun = @(x,y,t) -cos(k*pi*x).*sin(k*pi*y).*sin(k*pi*t);
    vfun = @(x,y,t) -sin(k*pi*x).*cos(k*pi*y).*sin(k*pi*t);
    ffun = @(x,y,t)  k*pi*(-1./(cfun(x,y))+2).*sin(k*pi*x).*sin(k*pi*y).*sin(k*pi*t);
%     ffun = @(x,y,t)  zeros(size(x));
    
    %% generate quadrature points
    global Pq CqM Vq
    
    Nq = 2*N+1;
    [rq sq wq] = Cubature2D(Nq); % integrate u*v*c
    Vq = Vandermonde2D(N,rq,sq)/V;
    xq = Vq*x; yq = Vq*y;
    
    %% construct the projection matrix for nodal basis Pq
    Pq = V*V'*Vq'*diag(wq);
    
    %% construct matrix Cq
    Cq = cfun(xq,yq);
    
    %% construct the projection matrix for cfun to degree M
    %     VMq = Vandermonde2D(M,rq,sq);
    %     CqM = VMq*VMq'*diag(wq)*Cq;
    
    d = zeros(Np,1);
    skk = 1;
    for ii = 0:N
        for jj = 0:N-ii
            if (ii+jj <= MM)
                d(skk) = 1;
            end
            skk = skk + 1;
        end
    end
    F = V*(diag(d)/V);
    CqM = Vq*F*Pq*Cq;
    
%     % interp
%     if MM>0
%         [rM sM] = Nodes2D(MM); [rM sM] = xytors(rM,sM);
%         VNM = Vandermonde2D(N,rM,sM)/V;
%         VMN = Vandermonde2D(MM,r,s)/Vandermonde2D(MM,rM,sM);
%         CqM = Vq*VMN*cfun(VNM*x,VNM*y);
%     end
    
    %% initial condition
    p = pfun(x,y,0);
    u = ufun(x,y,0);
    v = vfun(x,y,0);
    
    pw = p;
    uw = u;
    vw = v;
    
    time = 0;
    
    %% Runge-Kutta residual storage
    resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);
    resuw = zeros(Np,K); resvw = zeros(Np,K); respw = zeros(Np,K);
    
    %% compute time step size
    CN = (N+1)*(N+2)/2; % trace inequality constant
    CNh = max(CN*max(Fscale(:)));
    dt = CFL*2/CNh;
    wJq = diag(wq)*(Vq*J);
    
    %% outer time step loop
    tstep = 0;
    
    while (time<FinalTime)
        if(time+dt>FinalTime), dt = FinalTime-time; end
        
        for INTRK = 1:5
            
            timelocal = time + rk4c(INTRK)*dt;
            
            f = ffun(x,y,timelocal);
            
            [rhsp, rhsu, rhsv] = acousticsRHS2D_compare(p,u,v,f,CqM);
            
            % initiate and increment Runge-Kutta residuals
            resp = rk4a(INTRK)*resp + dt*rhsp;
            resu = rk4a(INTRK)*resu + dt*rhsu;
            resv = rk4a(INTRK)*resv + dt*rhsv;
            
            % update fields
            u = u+rk4b(INTRK)*resu;
            v = v+rk4b(INTRK)*resv;
            p = p+rk4b(INTRK)*resp;
            
            [rhsp, rhsu, rhsv] = acousticsRHS2D_compare(pw,uw,vw,f,Cq);
            respw = rk4a(INTRK)*respw + dt*rhsp;
            resuw = rk4a(INTRK)*resuw + dt*rhsu;
            resvw = rk4a(INTRK)*resvw + dt*rhsv;
            uw = uw + rk4b(INTRK)*resuw;
            vw = vw + rk4b(INTRK)*resvw;
            pw = pw + rk4b(INTRK)*respw;
        end
        
        % Increment time
        time = time+dt; tstep = tstep+1;
        if mod(tstep,10)==0
            fprintf('on timestep %d out of %d\n',tstep,ceil(FinalTime/dt))
        end
    end
    
    p_exact = pfun(xq,yq,FinalTime);
    
    p_quadrature = Vq * p;
    p_wadg = Vq * pw;
    error_l2(i) = sqrt(sum(sum(wJq.*(p_exact-p_quadrature).^2)));
    error_diff(i) = sqrt(sum(sum(wJq.*(p_wadg-p_quadrature).^2)));
    error_c(i) = sqrt(sum(sum(wJq.*(Cq-CqM).^2)));
end

%%

loglog(h,error_l2,'o--')
hold on
loglog(h,error_diff,'x--')
% loglog(h,error_c,'x--')
r = min(MM+3,N+1);
% r = 2*M+1;
% loglog(h,.1*h.^r,'k--')
% loglog(h,.25*h.^(MM+1),'k--')
% legend('Err','Diff','Error in c2')
legend('Err','Diff')
ratel2 = diff(log(error_l2))/log(.5);
ratediff = diff(log(error_diff))/log(.5);
title(sprintf('est L2 rate = %f, est diff rate = %f\n',mean(ratel2(end)),mean(ratediff(end))))
% figure
% vv = p_wadg-p_quadrature;
% vv = Cq-CqM;
% color_line3(xq,yq,vv,vv,'.')

%%


function [rhsp, rhsu, rhsv] = acousticsRHS2D_compare(p,u,v,f,c2)

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

global Pq Vq
rhsp = Pq*(c2.*(Vq*rhsp));
return;
end