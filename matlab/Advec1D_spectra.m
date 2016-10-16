function Advec1D_spectra

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u
Globals1D;
N = 4;

% Generate simple mesh
K1D = 8;
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

vmapP(1) = vmapM(end); % make periodic
vmapP(end) = vmapM(1);


M = inv(V*V');

rp = linspace(-1,1,100)';
Vp = Vandermonde1D(N,rp)/V;
xp =  Vp*x;


%% make CG projection 

R = zeros((N+1)*K,N*K);
offset = 0;
for e = 1:K
    r = (1:N+1) + (e-1)*(N+1);
    c = (1:N+1) + offset;
    R(r,c) = eye(N+1);
    offset = offset + N;
end
R(:,end) = []; 
R((N+1)*K,1) = 1; % make periodic

% P = R*((R'*kron(diag(J(1,:)),M)*R) \ (R'*kron(diag(J(1,:)),M)));
P = diag(1./sum(R,2))*R*R';
% return

%%

if 1 %nargout > 0
    tau = 1;
    A1 = zeros(Np*K);
    u = zeros(Np,K);
    for i = 1:Np*K
        u(i) = 1;
        rhs = AdvecRHS1D(u,tau);
        u(i) = 0;
        A1(:,i) = rhs(:);
    end
    lam = eig(A1);
    
    tau = 10;
    A = zeros(Np*K);
    u = zeros(Np,K);
    for i = 1:Np*K
        u(i) = 1;
        rhs = AdvecRHS1D(u,tau);
        u(i) = 0;
        A(:,i) = rhs(:);
    end
    
    Am = eye(size(A));
    C = (N*N*K);
    Ap = A/C;
    for i = 1:50
        Am = Am + Ap;
        Ap = (A/C)*Ap/(i+1);
    end
    
    norm(expm(A/C)- Am,'fro')

    keyboard
    [W D] = eig(A);
    d = diag(D);
    % remove +/- doubles
    [~,p] = sort(abs(lam),'descend');
    W = W(:,p);
    d = d(p);
    
    plot(lam,'o')
    hold on
    plot(d,'^')
    d2 = eig(P*A);
    plot(d2,'x')
    return
    
    %         ids = find(real(d) > -10 & real(d) < -1e-4); % converge to CG modes
    ids = find(real(d) < -10); % penalized modes
    d = d(ids);
    W = W(:,ids);
    
    i = 6;
    subplot(1,2,1)
    plot(real(d),imag(d),'o')
    hold on;plot(real(d(i)),imag(d(i)),'rx')
    subplot(1,2,2)
    ww = W(:,i);
    plot(xp,Vp*reshape(real(ww),N+1,K))
    hold on
    plot(xp,Vp*reshape(imag(ww),N+1,K))
    return
    
end

M = inv(V*V');

rp = linspace(-1,1,25)';
Vp = Vandermonde1D(N,rp)/V;
xp=  Vp*x;

[rq wq] = JacobiGQ(0,0,N+1);
Vq = Vandermonde1D(N,rq)/V;
xq = Vq*x;

% Set initial conditions
uex = @(x) sin(pi*x);
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
        rhsu = AdvecRHS1D(u,1);
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u + rk4b(INTRK)*resu;
        %         u(:) = Q*Q'*u(:);
    end;
    
    % Increment time
    time = time+dt;
    
    if mod(tstep,10)==0
        
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

err = diag(wq)*abs(Vq*u-uex(xq-time)).^2*diag(J(1,:));
err = sqrt(sum(err(:)))

P = Q*((Q'*kron(eye(K),M)*Q)\(Q'*kron(eye(K),M))); % L2 project onto Q
uQ = reshape(P*u(:),Np,K);
err = diag(wq)*abs(Vq*uQ-uex(xq-time)).^2*diag(J(1,:));
err = sqrt(sum(err(:)))
keyboard

function [rhsu] = AdvecRHS1D(u,alpha)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form field differences at faces
%alpha = 0;
du = zeros(Nfp*Nfaces,K);
du(:) = (u(vmapM)-u(vmapP)).*(nx(:)-alpha*abs(nx(:)))/2;

% compute right hand sides of the semi-discrete PDE
rhsu = -rx.*(Dr*u) + LIFT*(Fscale.*(du));
return

