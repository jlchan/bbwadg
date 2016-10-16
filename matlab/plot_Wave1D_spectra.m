% wave

function plot_Wave1D_spectra

% Driver script for solving the 1D advection equations
Globals1D;

% Order of polymomials used for approximation
N = 4;

% Generate simple mesh
K1D = 4;
[Nv, VX, K, EToV] = MeshGen1D(-1.0,1.0,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

% vmapP(1) = vmapM(end); % hack for periodic
% vmapP(end) = vmapM(1); % hack for periodic

rp = linspace(-1,1,500)';
Vp = Vandermonde1D(N,rp)/V;
xp = Vp*x;

%%
dhist = [];
dtau = 1e-1;
tauvec = 0:dtau:100;
tauvec = 100*exp(tauvec)/exp(100);
for ii = 1:length(tauvec)
    tau = tauvec(ii);
    
    Q = zeros(N+1,2*K);
    A = zeros(2*(N+1)*K);
    for i = 1:2*(N+1)*K
        Q(i) = 1;
        p = reshape(Q(:,1:K),N+1,K);
        u = reshape(Q(:,K+1:2*K),N+1,K);
        [rhsp rhsu] = WaveRHS1D(p,u,tau);
        Q(i) = 0;
        A(:,i) = [rhsp(:); rhsu(:)];
    end
    [W D] = eig(A); 
    d = diag(D);
    
    points{ii} = [real(d) imag(d)];
    
    if mod(ii,50) == 0
        disp(sprintf('on eig %d out of %d\n',ii,length(tauvec)))
    end
end

max_linking_distance = 3;
max_gap_closing = Inf;
debug = true;
disp('running simple tracker')
[ tracks, adjacency_tracks ] = simpletracker(points,...
    'MaxLinkingDistance', max_linking_distance, ...
    'MaxGapClosing', max_gap_closing, ...
    'Debug', debug);


plot(points{1}(:,1),points{1}(:,2),'o','markersize',9,'DisplayName','\tau = 0')
hold on
[~,id] = min(abs(tauvec-1));
plot(points{id}(:,1),points{id}(:,2),'^','markersize',9,'DisplayName','\tau = 1')
plot(points{end}(:,1),points{end}(:,2),'x','markersize',9,'DisplayName',sprintf('\\tau = %1.1f',tauvec(end)))

% exact eigs
lam_ex = pi/2*(1:Np*K); lam_ex = [lam_ex -lam_ex]; 
plot(zeros(size(lam_ex)),lam_ex,'s','markersize',11,'DisplayName','Exact eigenvalues')

legend show
n_tracks = numel(tracks);
colors = hsv(n_tracks);
all_points = vertcat(points{:});
for i_track = 1 : n_tracks
   
    % We use the adjacency tracks to retrieve the points coordinates. It
    % saves us a loop.
    
    track = adjacency_tracks{i_track};
    track_points = all_points(track, :);
    
    plot(track_points(:,1), track_points(:, 2), 'Color', colors(i_track, :),'linewidth',2)    
end
axis equal
set(gca,'fontsize',15)
print(gcf,'-dpng','trackedEigs.png')
keyboard
return



%% Solve Problem
% Set initial conditions
% p = sin(pi*x);
p = exp(-100*x.^2);
u = zeros(size(x));

FinalTime = 10;

time = 0;

% Runge-Kutta residual storage
resp = zeros(Np,K);
resu = zeros(Np,K);

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
dt  = .5*xmin;
% dt = min(dt,1/smax);
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

% tau = 1;

sk = 1;
% outer time step loop
for tstep=1:Nsteps
    for INTRK = 1:5
        [rhsp rhsu] = WaveRHS1D(p,u,tau);
        rhs = P*[rhsp(:); rhsu(:)];
        rhsp(:) = rhs((1:(N+1)*K));
        rhsu(:) = rhs((1:(N+1)*K) + (N+1)*K);
        
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        p = p + rk4b(INTRK)*resp;
        u = u + rk4b(INTRK)*resu;                
    end;
    % Increment time
    time = time+dt;
    if mod(tstep,10)==0
        plot(xp,Vp*p,'.-');
        axis([-1 1 -1 1])
        drawnow
    end
end;

function [rhsp rhsu] = WaveRHS1D(p,u,tau)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form field differences at faces
dp = zeros(Nfp*Nfaces,K);
du = zeros(Nfp*Nfaces,K);
dp(:) = p(vmapP)-p(vmapM);
du(:) = u(vmapP)-u(vmapM);

dp(mapB) = -2*p(vmapM(mapB));

pflux = .5*(tau*dp - du.*nx);
uflux = .5*(tau*du.*nx - dp).*nx;

% compute right hand sides of the semi-discrete PDE
rhsp = -rx.*(Dr*u) + LIFT*(Fscale.*pflux);
rhsu = -rx.*(Dr*p) + LIFT*(Fscale.*uflux);

return

