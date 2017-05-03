function [p x y] = save_solution_driver(Nin,K1D,cfun,FinalTime,c_flag)

% clear all, clear
clear -global *

Globals2D

if nargin==0
    Nin = 4;
    K1D = 16;
    c_flag = 0; % 0 = low storage
    a = 2;
    cfun = @(x,y) 1 + .5*sin(a*pi*x).*sin(a*pi*y) + (y>0);
%     cfun =@(x,y) 1 + .5*sin(2*pi*x + 1).*sin(pi*y + 2);
    FinalTime = .5; 
end
N = Nin;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

StartUp2D;

% keyboard
[rq sq w] = Cubature2D(2*N+1); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

% keyboard

%% params setup

global c cq W Vq Pq Minvc Mref xq yq USE_C_MASS

USE_C_MASS = c_flag;

% c = (V*V')*Vq'*(diag(w)*cfun(xq,yq)); % L2 projection - can have discontinuities
% cq = Vq*c;
cq = cfun(xq,yq);

% alias
Nuint = N-1;
[ruint suint] = Nodes2D(Nuint); [ruint suint] = xytors(ruint,suint);
V1 = Vandermonde2D(N,ruint,suint)/V;
V2 = Vandermonde2D(Nuint,r,s)/Vandermonde2D(Nuint,ruint,suint);
x1 = V1 * x;
y1 = V1 * y;
c = V2*cfun(x1,y1);

cq = Vq*c;
% keyboard

disp(sprintf('min continuous c value = %e\n',min(cq(:))))
if min(cq(:))<1e-12
    disp('negative velocity')
    keyboard
end


Pq = V*V'*Vq'*diag(w); % J's cancel out
Mref = inv(V*V');

k = 1; W = (2*k-1)/2*pi; p = cos(W*x).*cos(W*y); % frequency of solution
if nargin==0
    p = exp(-100*(x.^2 + (y-.25).^2));
end
u = zeros(Np, K); v = zeros(Np, K);

for e = 1:K
    Minvc{e} = Vq'*diag(w.*1./cq(:,e))*Vq;   
end



%%

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
dt = .5/(max(cq(:))*CN*max(Fscale(:)));

% outer time step loop
tstep = 0;
while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
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
    
    % Increment time
    time = time+dt; tstep = tstep+1;
    
    if mod(tstep,100)==0
        disp(sprintf('on time %f\n',time))
    end
    
    if nargin==0 && mod(tstep,5)==0
        clf                
        %         color_line3(xp,yp,vv,vv,'.');
        % surf(xp,yp,vv); shading interp; axis equal
        PlotField2D(N+1, x, y, p);        
        view(2)
        %         axis([-1 1 -1 1 -1.25 1.25]);
        title(sprintf('Time = %f\n',time))
        drawnow
    end
    
end
% keyboard
% return




function [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;
global c cq W Vq Pq Minvc Mref xq yq USE_C_MASS

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

% Impose reflective boundary conditions (p+ = -p-)
du(mapB) = 0; dv(mapB) = 0;
dp(mapB) = -2*p(vmapB);

% evaluate upwind fluxes
ndotdU = nx.*du + ny.*dv;
tau = 1;
fluxp =  tau*dp - ndotdU;
fluxu =  (tau*ndotdU - dp).*nx;
fluxv =  (tau*ndotdU - dp).*ny;

% local derivatives of fields
pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
rhsp =  -divU + LIFT*(Fscale.*fluxp)/2.0;

% % add forcing term
% f = sqrt(2)*W*(1-1./cq).*cos(W.*xq).*cos(W.*yq).*sin(sqrt(2)*W*time); % note that 2nd order form forcing is dfdt
% rhsp = rhsp + Pq*f;

if USE_C_MASS==1 % actual inversion of c-mass matrix
    for e = 1:K
        rhsp(:,e) = Minvc{e}\(Mref*rhsp(:,e));
    end
else  % projection
    rhsp = Pq*(cq.*(Vq*rhsp));
end

rhsu =  -dpdx + LIFT*(Fscale.*fluxu)/2.0;
rhsv =  -dpdy + LIFT*(Fscale.*fluxv)/2.0;

return;

