clear

Globals1D;

smoothKnots = 25;

% for K = [2 4 8 16];
sk = 1;
Nvec = 1:8;
for N = Nvec
    
    K = 2*N;%ceil(N/2);
    
    % Generate simple mesh
    [Nv, VX, K, EToV] = MeshGen1D(-1,1,K);
    
    if smoothKnots
        for i = 1:smoothKnots
            % squeeze control pts in to reduce CFL restriction
            re = linspace(-1,1,N+K)';
            reKsub = linspace(-1,1,K+1)';
            Ve = bspline_basismatrix(N+1, [VX(1)*ones(1,N) VX VX(end)*ones(1,N)], reKsub);
            VX = Ve*re; VX = VX(:)';
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
    
    t = [VX(1)*ones(1,N) VX VX(end)*ones(1,N)];
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
    
    %     if (N==4)
    %         rp = linspace(-1,1,500)';
    %         Vp = bspline_basismatrix(N+1,t,rp);
    %         keyboard
    %     end
    sk = sk + 1;
    if sk > 9
        sk
    end
end

C = max(CT/2,CM);
NN = 1:N;
figure
CTpoly = (NN+1).*(NN+1)/2;
plot(NN,CTpoly,'o--','linewidth',2,'markersize',8)
hold on
% plot(NN,CIpoly,'o--','linewidth',2,'markersize',8)
% plot(Nvec,C,'o--','linewidth',2,'markersize',8)
% plot(NN,1.5*NN,'k--')
plot(Nvec,CM,'s--','linewidth',2,'markersize',8)
plot(Nvec,CT,'x--','linewidth',2,'markersize',8)
legend('C_T poly','C_I poly','C_I','C_T')
axis on
xlabel('Degree N','fontsize',14)

set(gca,'fontsize',14)

% CT(3:end)
% CM(3:end)
% return
print_pgf_coordinates(NN,CTpoly)
print_pgf_coordinates(NN,CT)
print_pgf_coordinates(NN,CM)
