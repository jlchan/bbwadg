clear

Globals1D;

% N = 4;
K1D = 4;
sk = 1;
% Generate simple mesh
for N = 4
    [Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
    
    %     a = 0;
    %     VX = JacobiGL(0,0,K1D)';
    %     VX = VX*a + (1-2*acos(VX)/pi)*(1-a); % squeeze control pts in to reduce CFL restriction
    % VX(2:end-1) = (2/(N+1))*VX(2:end-1);
    
    % Initialize solver and construct grid and metric
    StartUp1D;
    
    vmapP(1) = vmapM(end); % make periodic
    vmapP(end) = vmapM(1);
    
    rp = linspace(-1,1,100)';
    [rq wq] = JacobiGQ(0,0,N);
        
    Vp = Vandermonde1D(N,rp)/V; xp = Vp*x;
        
    Vq = Vandermonde1D(N,rq)/V;
    Vrq = GradVandermonde1D(N,rq)/V;
    xq = Vq*x;
    wqJ = repmat(wq,1,K).*(Vq*J); %wqJ = wqJ(:);
    
    % wrt BB
    [VB] = bern_basis_1D(N,r);    
    [Vq Vrq] = bern_basis_1D(N,rq);
    Dr = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq)*Vrq);
    
    % BB mass/deriv matrices
    M = kron(diag(J(1,:)),Vq'*diag(wq)*Vq);
    S = kron(diag(J(1,:).*rx(1,:)),Vq'*diag(wq)*Vrq);
    %     Dx = M\S;
    Dx = kron(diag(rx(1,:)),Dr);
    
    t = [VX(1)*ones(1,N) VX VX(end)*ones(1,N)];
    % t = sort(cos(t*pi));
    B = bspline_basismatrix(N+1,t,x(:));
    
%     R = kron(eye(K),V)\B;
    R = kron(eye(K),VB)\B;
    
    dnormDG(sk) = max(abs(eig(Dx))); % block dmat matrix over [-1,1]: DG basis
    MB = (R'*M*R);
    dnormB(sk)= max(abs(eig(MB\(R'*S*R))));
    sk = sk + 1;
    
    plot(xp(:),kron(speye(K),Vp)*B)
    hold on;plot(t,ones(size(t)),'o')
end



% kron(speye(K),Vp)*B)

return
figure;
plot(dnormDG,'*--')
hold on
plot(dnormB,'o--')
% C = polyfit(1:N,dnormB,1);
% plot((1:N)*C(1),'rx--')

