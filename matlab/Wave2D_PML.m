function L2err = Wave2D_PML(Nin,K1D,c_flag,cfun)

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

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

StartUp2D;

[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

Nq = 2*N+1;
[rq sq w] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

%% params setup

global cq Vq Pq sigmax sigmay

cq = cfun(xq,yq);

disp(sprintf('min continuous c value = %e\n',min(cq(:))))
if min(cq(:))<1e-12
    disp('negative velocity')
    keyboard
end

Pq = V*V'*Vq'*diag(w); % J's cancel out

k = 1; % frequency of solution
W = (2*k-1)/2*pi;
p = cos(W*x).*cos(W*y);
p = exp(-100*(x.^2 + (y-0).^2));
u = zeros(Np, K); v = zeros(Np, K);

%% pml terms

showSol = 1;
delta = .5;
smax = 100;

M = 3;
L = 1-delta;
sigmax = (x > L).*(x-L).^M - (x < -L).*(x+L).^M;
sigmay = (y > L).*(y-L).^M - (y < -L).*(y+L).^M;

Q = zeros(Np,K);
R = zeros(Np,K);
resQ = zeros(Np,K); 
resR = zeros(Np,K);

P1 = eye(Np);
for e = 1:K
    %     sigmax(:,e) = max(sigmax(:,e));
    %     sigmay(:,e) = max(sigmay(:,e));
    sigmax(:,e) = mean(sigmax(:,e));
    sigmay(:,e) = mean(sigmay(:,e));
end
% P1 = V*diag([1;0;0;zeros(Np-3,1)])*inv(V); % reduce to const or linear modes
% 
sk = 1;
for e = 1:K
    if (max(sigmax(:,e)) > 1e-8 || max(sigmay(:,e))> 1e-8)
        pmlK(sk) = e;
        sk = sk + 1;
    end
end


% clf
% plot3(x,y,sigmax,'o'); return
% PlotMesh2D; axis on; return

sigmax = smax*sigmax/max(sigmax(:));
sigmay = smax*sigmay/max(sigmay(:));

M = inv(V*V');
%%
FinalTime = 25.0;

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
dt = 1/(max(cq(:))*CN*max(Fscale(:)));

dt = min(dt,1/(2*smax));

% outer time step loop
tstep = 1;
sk = 1;
% figure
while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        timelocal = time + rk4c(INTRK)*dt;
        
        [rhsp, rhsu, rhsv, rhsQ, rhsR] = acousticsRHS2D(p,u,v,Q,R);
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        resQ = rk4a(INTRK)*resQ + dt*rhsQ;
        resR = rk4a(INTRK)*resR + dt*rhsR;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        p = p+rk4b(INTRK)*resp;
        Q = Q+rk4b(INTRK)*resQ;
        R = R+rk4b(INTRK)*resR;
        
        p(:,pmlK) = P1*p(:,pmlK);
        u(:,pmlK) = P1*u(:,pmlK);
        v(:,pmlK) = P1*v(:,pmlK);
        Q = P1*Q;        
        R = P1*R;        
    end;
    
    if showSol && mod(tstep,10)==0
        clf
        %         vv = Vp*p;
        %         color_line3(x,y,p,p,'.');
        %         axis([-1 1 -1 1 -1.25 1.25]);
        PlotField2D(N+1, x, y, p); view(2)
        %         caxis([-.5 .5])
        colorbar
        drawnow
    elseif mod(tstep,10)==0
        Mp = M*p*diag(J(1,:));
        pK = p(:,pmlK); Mp = Mp(:,pmlK);
        l2norm(sk) = sqrt(pK(:)'*Mp(:));        
        tvec(sk) = time;
        semilogy(time,l2norm(sk),'.')
        hold on
        drawnow
        sk = sk + 1;
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
end

plot(tvec,l2norm);
return

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



function [rhsp, rhsu, rhsv rhsQ rhsR] = acousticsRHS2D(p,u,v,Q,R)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;
global cq Vq Pq sigmax sigmay

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

% Impose reflective boundary conditions (p+ = -p-)
% du(mapB) = 0; dv(mapB) = 0; dp(mapB) = -2*p(vmapB);
 du(mapB) = -2*u(vmapB); 
 dv(mapB) = -2*v(vmapB); 
 dp(mapB) = 0;

% evaluate upwind fluxes
ndotdU = nx.*du + ny.*dv;
tau = 1;
fluxp =  .5*(tau*dp - ndotdU);
fluxu =  .5*(tau*ndotdU - dp).*nx;
fluxv =  .5*(tau*ndotdU - dp).*ny;

pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
rhsp =  -divU + LIFT*(Fscale.*fluxp);
rhsu =  -dpdx + LIFT*(Fscale.*fluxu);
rhsv =  -dpdy + LIFT*(Fscale.*fluxv);

% pml terms
dQ = zeros(Nfp*Nfaces,K); dQ(:) = Q(vmapP)-Q(vmapM);
dR = zeros(Nfp*Nfaces,K); dR(:) = R(vmapP)-R(vmapM);
% dQ(mapB) = 0; 
% dR(mapB) = 0;

dQdy = ry.*(Dr*Q) + sy.*(Ds*Q);
dRdx = rx.*(Dr*R) + sx.*(Ds*R);
dRx = .5*(dR.*nx);
dQy = .5*(dQ.*ny);

% sxq = Vq*sigmax;
% syq = Vq*sigmay;
% pq = (Vq*p);
% uq = (Vq*u);
% vq = (Vq*v);
% rhsp = rhsp - Pq*((sxq + syq).*pq);
% rhsp = rhsp - Pq*(sxq.*(Vq*(dQdy + LIFT*(Fscale.*dQy)))) - Pq*(syq.*(Vq*(dRdx + LIFT*(Fscale.*dRx))));
% rhsu = rhsu - Pq*(sxq.*uq);
% rhsv = rhsv - Pq*(syq.*vq);

rhsp = rhsp - (sigmax + sigmay).*p;
rhsp = rhsp - sigmax.*(dQdy + LIFT*(Fscale.*dQy)) - sigmay.*(dRdx + LIFT*(Fscale.*dRx));
rhsu = rhsu - sigmax.*u;
rhsv = rhsv - sigmay.*v;


rhsQ = v;
rhsR = u;

% % non-const media
% rhsp = Pq*(cq.*(Vq*rhsp));

return;

