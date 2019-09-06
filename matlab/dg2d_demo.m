Globals2D

N = 7;
K1D = 4;
CFL = 2;
FinalTime = 2;

global tau
tau = 1;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

% aa = .25/K1D;
% VX = VX+aa*randn(size(VX));
% VY = VY+aa*randn(size(VY));

StartUp2D

Lx = max(VX)-min(VX);
Ly = max(VY)-min(VY);
BuildPeriodicMaps2D(Lx,Ly);

mapP = reshape(mapP,Nfp*Nfaces,K);

global rxJ sxJ ryJ syJ nxJ nyJ
rxJ = rx.*J; sxJ = sx.*J;
ryJ = ry.*J; syJ = sy.*J;
nxJ = nx.*sJ;
nyJ = ny.*sJ;

% set up plotting points
[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x;
yp = Vp*y;

%% runge kutta coefficients

uex = @(x) exp(-25*(x.^2+y.^2));
% uex = @(x) abs(x)<.5;
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
CN = (N+1)*(N+2)/2; % scaling for N
dt = CFL*h/CN; % set time-step
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

figure(1)
hold off
resu = zeros(size(u));
for i = 1:Nsteps
    
    for INTRK = 1:5        
        rhs = advecRHS(u);        
        resu = rk4a(INTRK)*resu + dt*rhs;
        u = u + rk4b(INTRK)*resu;
    end
        
    if mod(i,25)==0 || i==Nsteps               
        vv = Vp*u;
        clf
        color_line3(xp,yp,vv,vv,'.')        
        axis equal
        drawnow
    end
end

function rhs = advecRHS(u)

Globals2D

global rxJ sxJ ryJ syJ nxJ nyJ
global tau

% DG rhs evaluation
uf = u(Fmask,:);
flux = .5*(uf(mapP) + uf).*(nxJ+nyJ) - tau*.5*abs(nxJ+nyJ).*(uf(mapP)-uf); % average

dudxJ = rxJ.*(Dr*u) + sxJ.*(Ds*u);
dudyJ = ryJ.*(Dr*u) + syJ.*(Ds*u);

rhs = dudxJ + dudyJ + LIFT*((flux - uf.*(nxJ+nyJ)));
rhs = -rhs./J;

end








