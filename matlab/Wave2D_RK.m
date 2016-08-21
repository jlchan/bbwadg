function L2err = Wave2D_RK(Nin,K1D,c_flag,cfun)

% clear all, clear
clear -global *

Globals2D

if nargin==0
    %     Nin = 4; K1D = 16;
    Nin = 3; K1D = 8; c_flag = 0;
    cfun = @(x,y) ones(size(x));
    % cfun = @(x,y) 1 + .5*sin(pi*x).*sin(pi*y); % smooth velocity
%     cfun = @(x,y) (1 + .5*sin(2*pi*x).*sin(2*pi*y) + (y > 0)); % piecewise smooth velocity
end
N = Nin;

% filename = 'Grid/Other/block2.neu';
% filename = 'Grid/Maxwell2D/Maxwell05.neu';
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

StartUp2D;

% PlotMesh2D; return

[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

Nq = 2*N+1;
[rq sq w] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

%% params setup

global cq Vq Pq Minvc Mref USE_C_MASS

USE_C_MASS = c_flag;

cq = cfun(xq,yq);

disp(sprintf('min continuous c value = %e\n',min(cq(:))))
if min(cq(:))<1e-12
    disp('negative velocity')
    keyboard
end

Pq = V*V'*Vq'*diag(w); % J's cancel out
Mref = inv(V*V');

k = 1; % frequency of solution
W = (2*k-1)/2*pi;
p = cos(W*x).*cos(W*y);
% p = exp(-100*(x.^2 + (y-.25).^2));
u = zeros(Np, K); v = zeros(Np, K);


for e = 1:K
    Minvc{e} = Vq'*diag(w.*1./cq(:,e))*Vq;
    Pc{e} = Pq*diag(cq(:,e))*Vq;
    MinvMcM{e} = inv(V*V')*((Vq'*diag(w.*cq(:,e))*Vq)\inv(V*V'));
end

% test comparison between inverting M and doing a c-projection
if 0
    f = @(xq,yq) xq;
    M = blkdiag(Minvc{:});
    b = Vq'*diag(w)*f(xq,yq);
    u1 = reshape(M\b(:),Np,K);
    u2 = Pq*(cq.*(Vq*(V*V'*b))); % project onto const-c space
    norm(u1-u2,'fro')
end



%% check eigenvalues of DG matrix 
if 0 && nargin==0
    e = zeros(3*Np*K,1);
    A = zeros(3*Np*K);
    for i = 1:3*Np*K
        e(i) = 1;
        ids = 1:Np*K;
        p = reshape(e(ids),Np,K);
        u = reshape(e(ids + Np*K),Np,K);
        v = reshape(e(ids + 2*Np*K),Np,K);
        [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,0);
        A(:,i) = [rhsp(:);rhsu(:);rhsv(:)];
        e(i) = 0;
        if mod(i,100)==0
            disp(sprintf('on column %d out of %d\n',i,3*Np*K))
        end
    end
    lam = eig(A);
    hold on;
    plot(lam,'o')
    title(sprintf('Largest real part = %e',max(real(lam))))
    %     keyboard            
end


%%
FinalTime = 5.0;

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
dt = 1/(max(cq(:))*CN*max(Fscale(:)));

% outer time step loop
tstep = 0;
% figure
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
        %         vv = Vp*p;
        %         color_line3(x,y,p,p,'.');
        %         axis([-1 1 -1 1 -1.25 1.25]);
        PlotField2D(N+1, x, y, p); view(2)
        
        drawnow
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
end

% return

%%

%pex = @(x,y,time) cos(W*x).*cos(W*y).*cos(W*sqrt(2)*sqrt( )*time);
pex = @(x,y,time) cos(W*x).*cos(W*y).*cos(W*sqrt(2)*time);
pex = pex(xq,yq,FinalTime);
L2err2 = 0;
for e = 1:K
    diff = Vq*p(:,e) - pex(:,e);
    L2err2 = L2err2 + J(1,e)*sum(w.*diff.^2);
end

L2err = sqrt(sum(L2err2));
% if nargin==0
%     title(sprintf('L2 error at time %f = %e\n',FinalTime,L2err))
% end
% keyboard



function [rhsp, rhsu, rhsv] = acousticsRHS2D(p,u,v,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;
global cq Vq Pq Minvc Mref USE_C_MASS

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

pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
rhsp =  -divU + LIFT*(Fscale.*fluxp)/2.0;
rhsu =  -dpdx + LIFT*(Fscale.*fluxu)/2.0;
rhsv =  -dpdy + LIFT*(Fscale.*fluxv)/2.0;

if USE_C_MASS==1 % actual inversion of c2-weighted mass matrix
    if abs(time)<1e-8
        fprintf('full mass\n')
    end
    for e = 1:K
        rhsp(:,e) = Minvc{e}\(Mref*rhsp(:,e));
    end
elseif USE_C_MASS==0  % projection
    if abs(time)<1e-8
        fprintf('low storage projection\n')
    end
    rhsp = Pq*(cq.*(Vq*rhsp));
    
end


return;
