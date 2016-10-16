function maxeig = advec1D_cfl(Nin)

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u
if nargin==0
    Nin = 4;
end
Globals1D;

% Order of polymomials used for approximation
N = Nin;

% Generate simple mesh
K1D = 2;
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

vmapP(1) = vmapM(end); % make periodic
vmapP(end) = vmapM(1);


M = inv(V*V');

rp = linspace(-1,1,100)';
Vp = Vandermonde1D(N,rp)/V;
xp =  Vp*x;

if 1 %nargout > 0
    A = zeros(Np*K);
    u = zeros(Np,K);
    for i = 1:Np*K
        u(i) = 1;
        rhs = AdvecRHS1D(u);
        u(i) = 0;
        A(:,i) = rhs(:);
    end
    [W D] = eig(A);
    lam = diag(D);
    
    % remove +/- doubles
    [~,p] = sort(abs(lam),'descend');
    W = W(:,p);
    lam = lam(p);
    
    if 0
        for i = 1:size(W,2)
            w = real(W(:,i));
            subplot(1,2,1)
            w = reshape(w,Np,K);
            %             w = V\w;
            %             w = repmat(w(N+1,:),N+1,1);
            Nderivs = 0;
            for ii = 1:Nderivs
                w = rx.*(Dr*w);
            end
            vp = Vp*w;
            vp = vp/max(abs(vp(:)));
            plot(xp,vp,'o-')
            
            
            subplot(1,2,2)
            plot(lam,'o')
            hold on
            plot(real(lam(i)),imag(lam(i)),'rx')
            plot(-2*(N+1)*ones(10,1),linspace(min(imag(lam)),max(imag(lam)),10),'--')
            hold off
            pause
        end
    end
    
    %skip = floor(.5*sqrt(N)*K);
    %         skip = nnz(real(lam) < -2*(N+1));
    %     skip = nnz(imag(lam) < -2*(N+1));
    skip = nnz(abs(lam) > 2*(N+1)*K);
    
%     plot(lam,'o')
%     hold on
%     plot(lam(skip+1:end),'x')
%     axis equal
%     pause
    
    %     Q = null(W(:,1:skip)'*kron(eye(K),M)); % basis which is L2 orthog to bad modes
        Q = W;
%     Q = W(:,1:skip);
%     Q = eye(size(Q,1)) - Q*Q'; % directly project out high modes
    %     keyboard
    rb = [-1;1];
    Vbdr = Vandermonde1D(N,rb);
    [Ub,Sb,Vb] = svd(Vbdr);
    Nb = numel(rb);
    Vbb = Vb(:,1:Nb); Sbb = Sb(:,1:Nb);
    Vbi = Vb(:,Nb+1:end); % interior modal dofs
    
    % Vbb = pinv(Vbdr); % nodal version of vertex-decoupled
    
    if 1
        for i = 1:size(Q,2);
            %         q = real(Q(:,i));
            q = Q(:,i);
            q = reshape(q,Np,K);
            
            q = q/max(abs(q(:)));
            qp = Vp*q;

            clf
            subplot(1,2,1)
            %             plot(xp,qp,'-')
            hold on
            %         plot(x(i),real(q(i)),'o')
            plot([-1 1],[0 0],'k--')
            
            % orthogonalized interior vs edge components
            qi = Vbi*Vbi'*(V\q);
            qb = Vbb*Vbb'*(V\q);
            %         qb1 = Vbb(:,1)*Vbb(:,1)'*(V\q);
            %         qb2 = Vbb(:,2)*Vbb(:,2)'*(V\q);
            % keyboard
            qb1 = qb(:,1);
            qb2 = qb(:,2);
            plot(xp,Vandermonde1D(N,rp)*qi,'b')
            hold on
            plot(xp,Vandermonde1D(N,rp)*qb,'r')
            %         plot(xp,Vandermonde1D(N,rp)*qb1,'b')
            %         plot(xp,Vandermonde1D(N,rp)*qb2,'r')
            plot(x([1 N+1],:),Vandermonde1D(N,[-1; 1])*qb,'ks')
            
            %             plot(x(i),0,'*')
            
            subplot(1,2,2)
            imagesc(abs(Q'*kron(eye(K),M)*Q))
            hold on;plot(i,i,'ro','markersize',12)
            pause
        end
    end
    plot(real(lam),imag(lam),'o')
    hold on
    mu = eig(Q'*A*Q);
    plot(real(mu),imag(mu),'x')
    %     keyboard
    maxeig = max(abs(lam));
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
        rhsu = AdvecRHS1D(u);
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

