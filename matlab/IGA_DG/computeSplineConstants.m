% p-refinement

clear

Globals1D;

smoothKnots = 0;'opt';

% for K = [2 4 8 16];
sk = 1;
Nvec = 1:7; 12;
for N = Nvec
    
    K = ceil(N/2);
    K = N;
    K = 2*N;
    
    % Generate simple mesh
    [Nv, VX, K, EToV] = MeshGen1D(-1,1,K);
    t0 = [VX(1)*ones(1,N) VX VX(end)*ones(1,N)];
    t = t0;
        
    if strcmp(smoothKnots,'opt')
        if N==1
            VX = linspace(-1,1,K+1); % optimal for N=1
            t0 = [VX(1)*ones(1,N) VX VX(end)*ones(1,N)];
            t = t0;    
        else
            if K==2*N || K==N
                load optKnots.mat
                VX = optKnots{NBlist==N & KsubList==K};
                t0 = [VX(1)*ones(1,N) VX VX(end)*ones(1,N)];
                t = t0;                
            end
        end
    elseif smoothKnots > 0
        %for i = 1:smoothKnots
        tdiff = 1;
        iter = 1; maxit = 50;
        while tdiff > 1e-3 || iter < maxit
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
    
%     CIpoly(sk) = sqrt(max(eig(Vrq'*diag(wq)*Vrq,Vq'*diag(wq)*Vq)));
    
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
    CT(sk) = max(lam);
    
    lam = eig(R'*S*R,MB);
    CM(sk) = sqrt(max(lam));
    
    lam = eig(Vrq'*diag(wq)*Vrq,Vq'*diag(wq)*Vq);
    CMpoly(sk) = sqrt(max(lam));
    
    Ndofs(sk) = N+K;
    Kvec(sk) = K;
    sk = sk + 1;
    if sk > 1
        sk
    end
end

maxC = max(max(CMpoly))+5;

NN = Nvec;
CTpoly = (NN+1).*(NN+1)/2;
Ndofp = (NN+1);
% Ndofs = Ndofs;

figure
plot(Kvec,CTpoly,'o--')
hold on
plot(Kvec,CT./Kvec,'x--')

figure
plot(Kvec,CMpoly,'o--')
hold on
plot(Kvec,CM./Kvec,'x--')

return

% figure
figure(1)
plot(Nvec,CTpoly./Ndofp,'o--','linewidth',2,'markersize',8)
hold on
plot(Nvec,CT./Ndofs,'x--','linewidth',2,'markersize',8)
legend('C_T poly','C_T spline')
grid on
% axis([0 N+1 0 maxC])
xlabel('Degree N','fontsize',14)
set(gca,'fontsize',14)

figure(2)
plot(Nvec,CMpoly./Ndofp,'o--','linewidth',2,'markersize',8)
hold on
% plot(Nvec,C,'o--','linewidth',2,'markersize',8)
% plot(NN,1.5*NN,'k--')
plot(Nvec,CM./Ndofs,'x--','linewidth',2,'markersize',8)
legend('C_I poly','C_I spline')
grid on
% axis([0 N+1 0 maxC])
xlabel('Degree N','fontsize',14)
set(gca,'fontsize',14)

% CT(3:end)
% CM(3:end)
% return
%
disp('trace ineq')
print_pgf_coordinates(NN,CMpoly)
print_pgf_coordinates(NN,CM)

disp('inverse ineq')
print_pgf_coordinates(NN,CTpoly)
print_pgf_coordinates(NN,CT)

%% h-refinement

clear

Globals1D;

smoothKnots = 0;'opt';

% for K = [2 4 8 16];
sk = 1;
N = 4;
Kvec = 1:32;
for K = Kvec
        
    % Generate simple mesh
    [Nv, VX, K, EToV] = MeshGen1D(-1,1,K);
    t0 = [VX(1)*ones(1,N) VX VX(end)*ones(1,N)];
    t = t0;
        
    if strcmp(smoothKnots,'opt')
        if N==1
            VX = linspace(-1,1,K+1); % optimal for N=1
            t0 = [VX(1)*ones(1,N) VX VX(end)*ones(1,N)];
            t = t0;    
        else
            if K==2*N || K==N
                load optKnots.mat
                VX = optKnots{NBlist==N & KsubList==K};
                t0 = [VX(1)*ones(1,N) VX VX(end)*ones(1,N)];
                t = t0;                
            end
        end
    elseif smoothKnots > 0
        %for i = 1:smoothKnots
        tdiff = 1;
        iter = 1; maxit = 50;
        while tdiff > 1e-3 || iter < maxit
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
    
%     CIpoly(sk) = sqrt(max(eig(Vrq'*diag(wq)*Vrq,Vq'*diag(wq)*Vq)));
    
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
    CT(sk) = max(lam);%/K;
    
    lam = eig(R'*S*R,MB);
    CM(sk) = sqrt(max(lam));%/K;
    
    lam = eig(Vrq'*diag(wq)*Vrq,Vq'*diag(wq)*Vq);
    CMpoly(sk) = sqrt(max(lam));
    
    sk = sk + 1;
    if sk > 1
        sk
    end
end

plot(Kvec,CT,'o--')
hold on
plot(Kvec,CM,'x--')
