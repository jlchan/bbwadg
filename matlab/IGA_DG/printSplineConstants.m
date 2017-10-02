clear

Globals1D;

smoothKnots = 1;'opt';

% for K = [2 4 8 16];
sk = 1;
Nvec = 1:7;
for N = 2:4
        
    % Generate simple mesh
    [Nv, VX, K, EToV] = MeshGen1D(-1,1,K);
    t0 = [VX(1)*ones(1,N) VX VX(end)*ones(1,N)];
    t = t0;
        
    if strcmp(smoothKnots,'opt')
        load optKnots.mat
        t0 = optKnots{NBlist==N & KsubList==K};
        t = t0;        
    elseif smoothKnots > 0
        %for i = 1:smoothKnots
        tdiff = 1;
        iter = 1; maxit = 150;
        while tdiff > 1e-4 || iter < maxit
            % squeeze control pts in to reduce CFL restriction
            re = linspace(-1,1,N+K)';
            Ve = bspline_basismatrix(N+1, t, t0);
            tdiff = norm(t(:)-Ve*re);
            t = Ve*re;  
            iter = iter + 1;
        end
        if (iter > maxit)
            tdiff
        end
    end
    
    % Initialize solver and construct grid and metric
    StartUp1D;
    
    rp = linspace(-1,1,100)';
    [rq wq] = JacobiGQ(0,0,N);
    Vp = Vandermonde1D(N,rp)/V; xp = Vp*x;
    
    Vq = Vandermonde1D(N,rq)/V;
    Vrq = GradVandermonde1D(N,rq)/V;
    xq = Vq*x;
    wqJ = repmat(wq,1,K).*(Vq*J); %wqJ = wqJ(:);
    V = eye(Np);
        
    [VB] = bern_basis_1D(N,r);
    [Vq Vrq] = bern_basis_1D(N,rq);
    
    Dr = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq)*Vrq);
    
    % BB mass/deriv matrices
    M = kron(diag(J(1,:)),Vq'*diag(wq)*Vq);
    Sr = kron(eye(K),Vq'*diag(wq)*Vrq);
    S = kron(diag(rx(1,:)),Vrq'*diag(wq)*Vrq);
    
    B = bspline_basismatrix(N+1,t,x(:));
    
    R = kron(eye(K),VB)\B;
    R(abs(R)<1e-8) = 0;
    
    MB = (R'*M*R);
    
    e = zeros(size(MB,2));
    e(1,1) = 1;
    e(end,end) = 1;
    lam = eig(e,MB);
    CT(sk) = max(lam)/K;
    
    lam = eig(R'*S*R,MB);
    CM(sk) = sqrt(max(lam))/K;
    
    lam = eig(Vrq'*diag(wq)*Vrq,Vq'*diag(wq)*Vq);
    CMpoly(sk) = sqrt(max(lam));
    
    sk = sk + 1;
    if sk > 9
        sk
    end
end

maxC = max(max(CMpoly))+5;

NN = 1:N;
CTpoly = (NN+1).*(NN+1)/2;

% figure
figure(1)
plot(NN,CTpoly,'o--','linewidth',2,'markersize',8)
hold on
plot(Nvec,CT,'x--','linewidth',2,'markersize',8)
legend('C_T poly','C_T spline')
grid on
axis([0 N+1 0 maxC])
xlabel('Degree N','fontsize',14)
set(gca,'fontsize',14)


figure(2)
plot(NN,CMpoly,'o--','linewidth',2,'markersize',8)
hold on
% plot(Nvec,C,'o--','linewidth',2,'markersize',8)
% plot(NN,1.5*NN,'k--')
plot(Nvec,CM,'x--','linewidth',2,'markersize',8)
legend('C_I poly','C_I spline')
grid on
axis([0 N+1 0 maxC])
xlabel('Degree N','fontsize',14)
set(gca,'fontsize',14)

% CT(3:end)
% CM(3:end)
% return
print_pgf_coordinates(NN,CMpoly)
print_pgf_coordinates(NN,CTpoly)
print_pgf_coordinates(NN,CT)
print_pgf_coordinates(NN,CM)
