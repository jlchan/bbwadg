clear

K = 100;
N = 1; % polynomial degree
CFL = .75;
FinalTime = 2;

global Fmask mapP tau rxJ nxJ J M Dr L
tau = 0; % penalty parameter - try varying between [0,1]

%% set up reference element

% nodes + quadrature
r = JacobiGL(0,0,N); % interpolation points
[rq wq] = JacobiGL(0,0,N); % compute quadrature points

% VDM matrices
V = Vandermonde1D(N,r); % modal VDM
Vrq = GradVandermonde1D(N,rq)/V;
Vq = Vandermonde1D(N,rq)/V;

% mass and convection matrices
M = Vq'*diag(wq)*Vq;
Q = Vq'*diag(wq)*Vrq;

B = zeros(N+1,2);
B(1,1) = 1;
B(end,end) = 1;
Fmask = [1;N+1]; % face node indices

Dr = M\Q;
L = M\B;

%% physical elements

VX = linspace(-1,1,K+1);

x = zeros(N+1,K);
for e = 1:K
    h = VX(e+1)-VX(e);
    x(:,e) = h*(r+1)/2 + VX(e);
end

Nfp = 1;
Nfaces = 2;
mapP = reshape(1:Nfp*Nfaces*K,Nfp*Nfaces,K);
for e = 1:K
    mapP(:,e) = [(Nfp*Nfaces)*e-2, Nfp*Nfaces*(e+1)-1];
end
mapP(1) = 1; mapP(end) = Nfp*Nfaces*K;
mapP(end) = 1; mapP(1) = Nfp*Nfaces*K;
mapP = reshape(mapP,Nfp*Nfaces,K);

% geometric terms - these may be pedantic, but will set the stage for 2D/3D
nxJ = repmat([-1;1],1,K); % "scaled" normals
J = h/2*ones(N+1,K); % determinant of Jacobian of mapping
rx = 2/h*ones(N+1,K); % geometric mapping term for derivative
rxJ = rx.*J; % scaled geometric mapping

%% flux formulation

global fS
fS = @(uL,uR) (uL.^2 + uL.*uR + uR.^2)/6;

%% runge kutta coefficients

uex = @(x) exp(-25*x.^2);
% uex = @(x) abs(x)<.5;
uex = @(x) -sin(pi*x);
u = uex(x);

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

h = 2/K;
CN = (N+1)^2/2; % scaling for N
dt = CFL*h/CN; % set time-step
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

figure
resu = zeros(size(u));
for i = 1:Nsteps
    
    for INTRK = 1:5
        
        rhs = burgersFluxDiffRHS(u);
        
        if INTRK==5
            rhstest(i) = sum(sum(u.*J.*(M*rhs)));
        end        
        resu = rk4a(INTRK)*resu + dt*rhs;
        u = u + rk4b(INTRK)*resu;
    end
        
    if mod(i,10)==0 || i==Nsteps
        clf
        plot(x,u,'o--','linewidth',2)
        axis([-1 1 -1.5 1.5])
        title(sprintf('rhstest = %g, unorm = %g\n',rhstest(i),sum(sum(u.*(M*u)))))
        drawnow
    end
end

function rhs = burgersFluxDiffRHS(u)

global fS
global Fmask mapP tau rxJ nxJ J M Dr L

uf = u(Fmask,:);
uP = uf(mapP);

% entropy stable flux + Lax-Friedrichs penalty
Lf = abs(uf);
Lmax = max(Lf,Lf(mapP)); % "correct" way
% Lmax = Lf; % looks better but is wrong

flux = 1/6*(uP.^2 + uP.*uf + uf.^2).*nxJ - tau*.5*Lmax.*(uf(mapP)-uf);
rhs = L*((flux - uf.^2/2.*nxJ));

% volume contributions
for e = 1:size(u,2)
    [ux uy] = meshgrid(u(:,e));
    FS = fS(ux,uy);
    rhs(:,e) = rhs(:,e) + rxJ(:,e).*sum(2*Dr.*FS,2);
end 

% divide by Jacobian
rhs = -rhs./J;

end