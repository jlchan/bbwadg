clear

for Nin = 1:9
    
    Globals1D;
    
    % Order of polymomials used for approximation
    N = Nin;
    
    % Generate simple mesh
    K1D = N;
    [Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
    
    a = 0;
%     VX = JacobiGL(0,0,K1D)';
    VX = VX*a + (1-2*acos(VX)/pi)*(1-a); % squeeze control pts in to reduce CFL restriction
%     VX(2:end-1) = (2/(N+1))*VX(2:end-1);
    
    % Initialize solver and construct grid and metric
    StartUp1D;
    
    vmapP(1) = vmapM(end); % make periodic
    vmapP(end) = vmapM(1);
    
    rp = linspace(-1,1,100)';
    [rq wq] = JacobiGQ(0,0,N+4);
    
    Vp = Vandermonde1D(N,rp)/V; xp =  Vp*x;
    Vq = Vandermonde1D(N,rq)/V; xq = Vq*x;
    wqJ = repmat(wq,1,K).*(Vq*J); %wqJ = wqJ(:);
    
    [VB] = bern_basis_1D(N,r);
    [Vp Vpr] = bern_basis_1D(N,rp);
    [Vq Vqr] = bern_basis_1D(N,rq);
    Dr = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq)*Vqr);
    
    M = kron(diag(J(1,:)),Vq'*diag(wq)*Vq);
    Dx = kron(diag(rx(1,:)),Dr);
    S = M*Dx;
    
    
    t = [VX(1)*ones(1,N) VX VX(end)*ones(1,N)];
    % t = sort(cos(t*pi));
    B = bspline_basismatrix(N+1,t,x(:));
    R = kron(eye(K),VB)\B;
    
    dnormDG = norm(M\S); % block dmat matrix over [-1,1]: DG basis
    dnormB = norm((R'*M*R)\(R'*S*R));
    % h = VX(2)-VX(1);
    % .25*N^2/h
    
%     if (N > 5),plot(xp(:),kron(eye(K),Vp)*R,'o-'); hold on; plot(VX,VX*0,'s');return,end
    
    % test errors
    % clf
    f = @(x) exp(x);
    
    u = (R'*M*R) \ R'*(kron(eye(K),Vq))'*(wqJ(:).*f(xq(:)));
    u = reshape(R*u,Np,K);
    e = wqJ.*(f(xq)-Vq*u).^2;
    eB = sqrt(sum(e(:)));
    % plot(xp,abs(f(xp)-Vp*u))
    % hold on
    
    u = (Vq'*diag(wq)*Vq) \ (Vq'*diag(wq)*f(xq));
    e = wqJ.*(f(xq)-Vq*u).^2;
    eDG = sqrt(sum(e(:)));
        
    %     plot(xp,f(xp))
    %     hold on
    %     plot(xp,Vp*u,'.')
    %     % plot(xp,abs(f(xp)-Vp*u),'--')
    %     disp(sprintf('dnormDG = %f, dnormB = %f\n',dnormDG,dnormB))
    %     disp(sprintf('L2err DG = %g, L2err B = %g\n',eDG,eB))
    
    errDG(Nin) = eDG;
    errB(Nin) = eB;
    dnDG(Nin) = dnormDG;
    dnB(Nin) = dnormB;
end

semilogy(errDG)
hold on
semilogy(errB)

figure
clf
plot(dnDG)
hold on
plot(dnB)

