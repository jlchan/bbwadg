clear
Nvec = 2:5;
K1D = 1;
Kvec = 32:32:256;
smoothKnots = 75;
sk1 =1;
for NB = Nvec    
    sk2 = 1;
    for Ksub = Kvec
        rho(sk1,sk2) = IPDG_IGA(NB,Ksub,K1D,smoothKnots);
        sk2 = sk2 + 1;
    end
    
%     plot(Kvec,sqrt(rho),'o--')
%     hold on
%     h = 1./Kvec(:);
    C = [Kvec(:).^0 Kvec(:)]\(sqrt(rho(sk1,:)'));
    slope(NB) = C(2);
%     plot(Kvec,Kvec*C(2),'x--')
    sk1 = sk1+1;
end



function rho = IPDG_IGA(NBin,Ksub,K1D,smoothKnots)

% Driver script for solving IPDG wave equation
Globals1D;

global NB

% Order of polymomials used for approximation
if nargin==0
    NB = 4;
    Ksub = 32;
    K1D = 1;
    smoothKnots = 50;
else
    NB=NBin;
end


N = NB+Ksub-1;

FinalTime = 1;

dofs = (N+1)*K1D

[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

%% cubature test

[rq wq] = JacobiGQ(0,0,N+4);
Vq = Vandermonde1D(N, rq)/V;
xq = Vq*x;
wJq = diag(wq)*(Vq*J);

%% penalty constant stuff

global tau
if Ksub==1
    CN = (N+1)*(N+1); % 1D trace const
    tau = 2*CN*Fscale;  % trace const * sum(face areas)*(ref elem surface area)
else
    tau = 2*NB*Ksub*Fscale; % estimate it based on Ksub
end

%%
% for plotting - build coordinates of all the nodes
rp = linspace(-1,1,150);
Vp = Vandermonde1D(N,rp)/V;
xp = Vp*x;

%% set init condition

%uex = @(x,t) cos(pi*x/2);
uex = @(x,t) exp(-5^2*x.^2);
u = uex(x,0);

%% make splines
if 1
    [BVDM M Dr R rBq wBq Bq Brq] = bsplineVDM(NB,Ksub,r,smoothKnots); % VDM for interp, mass, M\S
    
    %     [Bq] = bsplineVDM(NB,Ksub,rq,smoothKnots); % VDM for interp, mass, M\S
    %     M = Bq'*diag(wq)*Bq; % under-integrate mass matrix
    
    Mf = zeros(N+1,2);
    Mf(1,1) = 1;
    Mf(N+1,2) = 1;
    LIFT = M\Mf;
    
    xq = Vandermonde1D(N,rBq)/V * x;
    rq = rBq;
    wq = wBq;
    wJq = diag(wq)*(Vandermonde1D(N,rBq)/V*J);
    Vq = Bq;
    
    Vp = bsplineVDM(NB,Ksub,rp,smoothKnots);
    
    % assume uniform mesh for simplicity
    rx = rx(1);
    
    % Set initial conditions
    Pq = M\(Bq'*diag(wBq));
    
    xBq = Vandermonde1D(N,rBq)/V * x;
    u = Pq*uex(xBq);
    %     u = BVDM\uex(x);
end

%% check matrix symmetry
if 1
    u = zeros(Np,K);
    A = zeros(Np*K);
    for i = 1:Np*K
        u(i) = 1;
        [rhsu] = IPDGRHS(u,0);
        u(i) = 0;
        
        A(:,i) = rhsu(:);
    end
    %     MM = kron(diag(J(1,:)),M); % global mass matrix
    %     A = MM*A; % remove implicit mass matrix inverse
    A(abs(A)<1e-8) = 0;
    
    %     A = Brq'*diag(wq)*Brq;
    %     A = -MM\A;
    %     % A(1,:) = 0; A(:,1) = 0; A(1,1) = 1;
    %     % A(end,:) = 0; A(:,end) = 0; A(end,end) = 1;
    %     A = A(2:end-1,2:end-1);
    
    rho = max(abs(eig(A)));
    return
end

%%

dt = 2./max(tau(:)); % empirical

v = zeros(Np,K);

% outer time step loop
time = 0;
tstep = 0;

totalSteps = floor(FinalTime/dt);

resu = zeros(Np,K);
resv = zeros(Np,K);

while (time<FinalTime)
    
    if(time+dt>FinalTime), dt = (FinalTime-time);  end;
    
    % low storage RK
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        rhsv = IPDGRHS(u,timelocal);
        rhsu = v;
        
        resv = rk4a(INTRK)*resv + dt*rhsv;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        
        v = v + rk4b(INTRK)*resv;
        u = u + rk4b(INTRK)*resu;
    end;
    
    time = time+dt;
    tstep = tstep+1; % Increment time
    
    if (mod(tstep,1)==0) || abs(time-FinalTime)<1e-8
        clf
        plot(xp,Vp*u)
        axis([-1 1 -1 1])
        title(sprintf('time t = %f',time))
        
        drawnow
    end
    
end

% keyboard
end

function [rhsu] = IPDGRHS(u,time)

Globals1D;

% Define field differences at faces
du  = zeros(Nfp*Nfaces,K); du(:) = u(vmapM)-u(vmapP);

% impose boundary condition - Dirichlet BC's
du(1) = 2*u(1); du(end) = 2*u(end); % reflection

% Compute q
fluxu = nx.*du/2.0;
ux = rx.*(Dr*u);
q = ux - LIFT*(Fscale.*fluxu);
dq = zeros(Nfp*Nfaces,K); dq(:) = q(vmapM)-(ux(vmapM)+ux(vmapP))/2.0;

% impose boundary condition - Neumann BC's
dq(1) = q(1) - ux(1);
dq(end) = q(end) - ux(end);

% evaluate fluxes
global NB
hmin = 2.0/max(max(rx)); tau = (NB+1)^2/hmin;
fluxq = nx.*(dq + tau*nx.*du);

% compute right hand sides of the semi-discrete PDE
rhsu = rx.*(Dr*q) - LIFT*(Fscale.*fluxq);
% rhsu(1) = 0;
% rhsu(end) = 0;
end
