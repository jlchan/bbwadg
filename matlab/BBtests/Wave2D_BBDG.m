function Wave2D


Globals2D

N = 7;
K1D = 8;
c_flag = 0;
FinalTime = .5;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
StartUp2D;

% BuildPeriodicMaps2D(1,1);

% PlotMesh2D; return

[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

%% init cond setup

x0 = 0; y0 = .1;
p = exp(-200*((x-x0).^2 + (y-y0).^2));
u = zeros(Np, K); 
v = zeros(Np, K);

%% using BB basis

if 1 % set to 1 to use BB
    
    [V Vr Vs V1 V2 V3] = bern_basis_tri(N,r,s);
    Dr = V\Vr;
    Ds = V\Vs;
    
    % plotting interp matrix
    Vp = bern_basis_tri(N,rp,sp);
    
    % map nodal coeffs to BB coeffs
    p = V\p;
    u = V\u;
    v = V\v;
    
    % transform lift
    r1D = JacobiGL(0,0,N);
    V1D = bern_basis_1D(N,r1D); % assume GLL nodes on edge
    LIFT = V\(LIFT*blkdiag(V1D,V1D,V1D));
    
    % ======= alternative construction of LIFT        
    
    % degree elev from deg N to N+1
    ENp1 = bern_basis_1D(N+1,JacobiGL(0,0,N+1)) \ bern_basis_1D(N,JacobiGL(0,0,N+1));  
    
    % lift reduction
    EL1 = [];
    for j = 0:N;
        cj = (-1)^j * nchoosek(N,j)/(1+j);
        
        Ej = bern_basis_1D(N,r1D)\bern_basis_1D(N-j,r1D);
        EL1 = [EL1; cj*Ej'];
    end
    L0 = ENp1'*ENp1 * (N+1)^2/2;
    
    % compute block cols of LIFT by symmetry of simplex
    pids = cell(Nfaces,1); % permutation ids
    ids = 1:Nfp;
    for f = 1:Nfaces
        for i = 1:Np
            Lrow = LIFT(i,ids + (f-1)*Nfp);
            diff = LIFT(:,ids)-repmat(Lrow,Np,1);
            [~,rowid] = min(sum(abs(diff),2));
            pids{f}(i) = rowid;
        end
    end
    EL = [EL1 EL1(pids{2},:) EL1(pids{3},:)];
    LIFT = EL * blkdiag(L0,L0,L0);    
    
end

%%

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
CNh = max(CN*max(Fscale(:)));
dt = 2/CNh;

% outer time step loop
tstep = 0;
figure
% colormap(gray)
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
    
    if 1 && nargin==0 && mod(tstep,10)==0
        clf
        pp = p;
        vv = Vp*pp;
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f',time))
        drawnow
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
end

% axis off
% view(0,0)


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

% Impose reflective boundary conditions (p+ = -p-)
ndotdU(mapB) = 0;
dp(mapB) = -2*p(vmapB);

% % free surface BCs
% ndotdU(mapB) = -2*(nx(mapB).*u(vmapM(mapB)) + ny(mapB).*v(vmapM(mapB)));
% dp(mapB) = 0;

% % basic ABCs
% fluxp(mapB) = (nx(mapB).*u(vmapM(mapB)) + ny(mapB).*v(vmapM(mapB)));
% fluxu(mapB) = p(vmapM(mapB)).*nx(mapB);
% fluxv(mapB) = p(vmapM(mapB)).*ny(mapB);

tau = 1;
fluxp =  tau*dp - ndotdU;
fluxu =  (tau*ndotdU - dp).*nx;
fluxv =  (tau*ndotdU - dp).*ny;

pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
rhsp =  -divU + LIFT*(Fscale.*fluxp)/2.0;
rhsu =  -dpdx + LIFT*(Fscale.*fluxu)/2.0;
rhsv =  -dpdy + LIFT*(Fscale.*fluxv)/2.0;




return;

