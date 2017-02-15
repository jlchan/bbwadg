clear

Globals1D;

% for K1D = [2 4 8 16];
sk = 1;
Nvec = 1:6;
for N = Nvec
    K1D = N;
    % Generate simple mesh
    [Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
    
    %     a = 0;
    %     VX = JacobiGL(0,0,K1D)';
    %     VX = VX*a + (1-2*acos(VX)/pi)*(1-a); % squeeze control pts in to reduce CFL restriction
    % VX(2:end-1) = (2/(N+1))*VX(2:end-1);
    
    % Initialize solver and construct grid and metric
    StartUp1D;
    
    Dx = kron(diag(rx(1,:)),Dr);
    
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
    S = kron(diag(rx(1,:)),Vrq'*diag(wq)*Vrq);
    %     Dx = M\S;
    
    t = [VX(1)*ones(1,N) VX VX(end)*ones(1,N)];
    % t = sort(cos(t*pi));
    B = bspline_basismatrix(N+1,t,x(:));
    
    R = kron(eye(K),VB)\B;
    
    MB = (R'*M*R);
    if 0
        f = @(x) sin(pi*x);
        
        b = Vq'*(wqJ.*f(xq));
        u = R*(MB\(R'*b(:)));
        % u = R*(B\f(x(:)));
        u = reshape(u,N+1,K);
        % plot(xp,Vp*VB*u,'--','linewidth',2)
        % hold on
        % plot(xp,f(xp),'-','linewidth',2)
        % hold on;plot(t,ones(size(t)),'o')
        
        err = wqJ.*(f(xq)-Vq*u).^2;
        L2err = sqrt(sum(err(:)));
        L2err * size(R,2)
        % kron(speye(K),Vp)*B)
    end
    e = zeros(size(MB,2));
    e(1,1) = 1;    
    e(end,end) = 1;
    lam = eig(e,MB);
    CT(sk) = max(lam)/K1D;
    
    lam = eig(R'*S*R,MB);
    CM(sk) = sqrt(max(lam))/K1D;
    sk = sk + 1
end

C = max(CT,CM);
plot(Nvec,C,'o--')
hold on
NN = 1:N;
plot(NN,(NN+1).*(NN+2)/2,'x--')
