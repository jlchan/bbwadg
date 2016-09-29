function advec1D_pp

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;

% Order of polymomials used for approximation
N = 3;

% Generate simple mesh
K1D = 64;
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

vmapP(1) = vmapM(end); % make periodic
vmapP(end) = vmapM(1);

rp = linspace(-1,1,25)';
Vp = Vandermonde1D(N,rp)/V;
xp=  Vp*x;

global VB DB VBe 

re = linspace(-1,1,N+1)';
[VB VBr] = bern_basis_1D(N,r);
DB = VB\VBr;
DB(abs(DB)<1e-8) = 0;
VBe = bern_basis_1D(N,re);
xe = Vandermonde1D(N,re)/V * x;

% Set initial conditions
uex = @(x) sin(2*pi*x);
u = uex(x);

% Solve Problem
FinalTime = 1;

time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K);
xmin = min(abs(x(1,:)-x(2,:))); % compute time step size
dt   = .75*xmin;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;

% outer time step loop
for tstep=1:Nsteps
    
    % low storage RK
    for INTRK = 1:5        
        rhsu = AdvecRHS1D(u);        
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u + rk4b(INTRK)*resu;        
    end;

    % Increment time
    time = time+dt;
    
    if 0 %mod(tstep,10)==0
        
        plot(xp,Vp*u,'-')
        hold on;
        plot(x,u,'o')
        plot(xp,uex(xp-time),'--');        
        hold off
        axis([-1 1 -1 3])
        title(sprintf('Time = %f\n',time))
        drawnow
    end
    
end;

% plot(xp,Vp*u,'-')
% hold on;
% plot(x,u,'o')
clf
eDG = Vp*u-uex(xp-time);

NB = 2*N+1;
rB = JacobiGL(0,0,NB);
xB = Vandermonde1D(N,rB)/V*x;
[VB] = bern_basis_1D(NB,rB);

t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % repeat knots at end
B = bspline_basismatrix(NB+1,t,xB(:));
% plot(xB(:),B); return
% plot(xB,0*xB,'o'); return
R = kron(eye(K),VB)\B;
R(abs(R)<1e-8) = 0;
% Vp = bern_basis_1D(NB,rp);clf;hold on
% for i = 1:size(R,2)
%     plot(xp,Vp*reshape(R(:,i),NB+1,K));    
% end
% return

[rBq wBq] = JacobiGQ(0,0,NB);
NBq = length(wBq);
wBqJ = repmat(wBq(:),1,K).*repmat(J(1,:),length(wBq(:)),1);
xBq = Vandermonde1D(N,rBq)/V*x;
VBq = bern_basis_1D(NB,rBq);

L2eDG = wBqJ.*(Vandermonde1D(N,rBq)/V*u-uex(xBq-time)).^2;
L2eDG = L2eDG(:,round(K/2)-2:round(K/2)+2); % look at middle err
L2eDG = sqrt(sum(L2eDG(:)));

uq = Vandermonde1D(N,rBq)/V * u; 

% uex = @(x) (x+time).^N; uq = uex(xBq-time); % test recover poly

Bq = kron(eye(K),VBq)*R;  % spline VDM
MB = Bq'*diag(wBqJ(:))*Bq; % spline mass matrix

% project solution onto global splines
upp = (MB\(Bq'*(wBqJ(:).*uq(:)))); 
upp = reshape(R*upp,NB+1,K); 

eDGP = bern_basis_1D(NB,rp)*upp-uex(xp-time); % projection error

L2eDGP = wBqJ.*(bern_basis_1D(NB,rBq)*upp-uex(xBq-time)).^2;
L2eDGP = L2eDGP(:,round(K/2)-2:round(K/2)+2); % look at middle err
L2eDGP = sqrt(sum(L2eDGP(:)));

% local projection onto BB basis
% MBlocal = VBq'*diag(wBq)*VBq;
% upp = (MBloca)\ (VBq'*diag(wBq)*reshape(uq,NBq,K)); 
% upp = pinv(R)*upp(:);%diag(1./sum(R',2))*R'*upp(:);
% norm(inv(MB) - pinv(R)*kron(eye(K),inv(MBloca))*R,'fro')
% upp = reshape(R*upp,NB+1,K); 

% bezier projection 
ids = (1:NBq);
for e = 1:K
    for i = 1:size(Bq,2)
        w(i,e) = Bq(ids + (e-1)*NBq,i)'*wBqJ(:,e); % get weights
    end
end
wdenom = Bq'*wBqJ(:); % integrals of all Bq
WW = diag(1./wdenom)*w; % WW(i,e) = int_e(Bi) / int_[-1,1] Bi
for e = 1:K    
    W(:,e) = WW(ids+(e-1),e);
end
for e = 1:K    
    r = (1:NBq) + NBq*(e-1);
    c = (1:NBq) + (e-1);        
    Blocal{e} = Bq(r,c);
end
Blocal = blkdiag(Blocal{:});
RB = Blocal\Bq; 
RB(abs(RB)<1e-8) = 0;

RBP = B*RB'*diag(W(:))/Blocal;

% mesh(log(abs(RBP*kron(eye(K),Vandermonde1D(N,rBq)/V )))) % shows a very localized smoother
upp = RBP*uq(:); % represet in spline basis
% keyboard
% upp(:) = Blocal \ uq(:); % representation in local spline pieces    
% upp = W.*upp; % apply weighting to local spline coeffs
% upp = RB'*upp(:); % assemble global splines from pieces
% upp = B*upp;
upp = bern_basis_1D(NB,rB)\reshape(upp,NBq,K); % represent in BB basis for postprocessing

eDGBP = bern_basis_1D(NB,rp)*upp-uex(xp-time); % bezier projection
L2eDGBP = wBqJ.*(bern_basis_1D(NB,rBq)*upp-uex(xBq-time)).^2;
L2eDGBP = L2eDGBP(:,round(K/2)-2:round(K/2)+2); % look at middle err
L2eDGBP = sqrt(sum(L2eDGBP(:)));

plot(xp,abs(eDG),'--');
hold on
plot(xp,abs(eDGP),'s-');
plot(xp,abs(eDGBP),'.');

title(sprintf('DG err = %g, proj DG err = %g, Bezier proj DG err = %g\n',L2eDG,L2eDGP,L2eDGBP))

return

semilogy(xp,abs(eDG),'-');
hold on
plot(xp,abs(eDGP),'--');
plot(xp,abs(eDGBP),'.--');

keyboard


function [rhsu] = AdvecRHS1D(u)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form field differences at faces
alpha = 0;
du = zeros(Nfp*Nfaces,K);
du(:) = (u(vmapM)-u(vmapP)).*(nx(:)-(1-alpha)*abs(nx(:)))/2;

% compute right hand sides of the semi-discrete PDE
rhsu = -rx.*(Dr*u) + LIFT*(Fscale.*(du));
return

