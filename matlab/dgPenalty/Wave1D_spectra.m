% wave

function Wave1D_spectra

% Driver script for solving the 1D advection equations
Globals1D;

% Order of polymomials used for approximation
N = 4;

% Generate simple mesh
K1D = 8;
[Nv, VX, K, EToV] = MeshGen1D(-1.0,1.0,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

vmapP(1) = vmapM(end); % hack for periodic
vmapP(end) = vmapM(1); % hack for periodic

rp = linspace(-1,1,500)';
Vp = Vandermonde1D(N,rp)/V;
xp = Vp*x;

%% make CG projection 

M = inv(V*V');
R = zeros((N+1)*K,N*K);
offset = 0;
for e = 1:K
    r = (1:N+1) + (e-1)*(N+1);
    c = (1:N+1) + offset;
    R(r,c) = eye(N+1);
    offset = offset + N;
end
R(:,end) = []; R((N+1)*K,1) = 1; % make periodic

P = R*((R'*kron(diag(J(1,:)),M)*R) \ (R'*kron(diag(J(1,:)),M)));
% P = eye((N+1)*K);
P = kron(eye(2),P);

% P = diag(1./sum(R,2))*R*R';
% return

%%
for tau = 0:.01:.5
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
    plot(d,'o')
    hold on
%     plot(tau,tau,'*')
    drawnow
end
return
% [~,p] = sort(real(d),'ascend');
% d = d(p);
% W = W(:,p);

% ids = find(real(d) > -1250 & real(d) < 1e-4);
% ids = find(real(d) < -50);
% d = d(ids);
% W = W(:,ids);


plot(real(d),imag(d),'o')
hold on

[W2 D2] = eig(P*A);
d2 = diag(D2);
plot(real(d2),imag(d2),'s')

% clf
% w = W2(:,d2 < -10); w = w(1:(N+1)*K);
% plot(xp,Vp*reshape(real(w),N+1,K))


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

% dp(mapB) = -2*p(vmapM(mapB));

pflux = .5*(tau*dp - du.*nx);
uflux = .5*(tau*du.*nx - dp).*nx;

% compute right hand sides of the semi-discrete PDE
rhsp = -rx.*(Dr*u) + LIFT*(Fscale.*pflux);
rhsu = -rx.*(Dr*p) + LIFT*(Fscale.*uflux);

return

