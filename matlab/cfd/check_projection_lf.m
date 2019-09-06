% Order of polymomials used for approximation
clear -globals
clear

Globals1D

N = 1;
sk = 1;
for i = 1:10
    K1D = 2^i;
    h(sk) = 2/K1D;
    
    [Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
    StartUp1D;
    
    [rq wq] = JacobiGL(0,0,N);
    Vq = Vandermonde1D(N,rq)/V;
    Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));
    
    [rq2 wq2] = JacobiGQ(0,0,N);
    Vq2 = Vandermonde1D(N,rq2)/V;
    Pq2 = (Vq2'*diag(wq2)*Vq2)\(Vq2'*diag(wq2));
    
    rp = linspace(-1,1,100)';
    Vp = Vandermonde1D(N,rp)/V;
    
    xq = Vq*x;
    xq2 = Vq2*x;
    xp = Vp*x;
    
    f = @(x) x > .1;
    f = @(x) (x > .1).*x;% + exp(x).*sin(pi*x)
    f = @(x) tanh(200*x);  
%     f = @(x) exp(1+sin(6*pi*x));
    
    err(sk) = sum(J(1,:).*(wq2'*(Vq2*(Pq*f(xq) - Pq2*f(xq2))).^2));
    
%     subplot(1,2,1)
%     plot(xp,Vp*Pq*f(xq))
%     hold on
%     plot(xp,Vp*Pq2*f(xq2))

%     plot(xp,Vp*(Pq*f(xq)-Pq2*f(xq2)))
%     hold on
%     title(sprintf('diff = %g\n',err(sk)))
    sk = sk + 1;
    
    % subplot(1,2,2)
    % bar([abs(Pq*f(rq)),abs(Pq2*f(rq2))],'grouped')
end

figure
loglog(h,err,'o--')
hold on
loglog(h,err(1)/h(1)*h,'k-')
loglog(h,err(1)/h(1)*h.^(2*N+2),'k-')