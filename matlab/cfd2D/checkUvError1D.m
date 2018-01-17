clear
Globals1D;

Kvec = [4 8 16 32 64];
h = 1./Kvec;
for N = 1:5;    
sk = 1;
for K1D = Kvec
    projectV = 1;
    
    
    useSBP = 0;
    
    FinalTime = 1.8;
    FinalTime = .2;
    FinalTime = 2;
    opt = 1;
    
    global tau
    tau = 0;
    
    r = JacobiGL(0,0,N);
    % r = JacobiGQ(0,0,N);
    
    [rq wq] = JacobiGL(0,0,N); rq = rq*(1-1e-11);
    % [rq wq] = JacobiGL(0,0,N+4); rq = rq*(1-1e-11);
    [rq wq] = JacobiGQ(0,0,N+1);
    
    Nq = length(rq);
    
    if opt==1 || opt==0
        [Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
    elseif opt==2
        [Nv, VX, K, EToV] = MeshGen1D(-.5,.5,K1D);
    elseif opt==3
        [Nv, VX, K, EToV] = MeshGen1D(-5,5,K1D);
    elseif opt==4
        [Nv, VX, K, EToV] = MeshGen1D(0,1,K1D);
    end
    
    % Initialize solver and construct grid and metric
    StartUp1D;
    
    V = Vandermonde1D(N,r);
    Dr = GradVandermonde1D(N,r)/V;
    Vq = Vandermonde1D(N,rq)/V;
    Vf = Vandermonde1D(N,[-1 1])/V;
    
    [rq2 wq2] = JacobiGQ(0,0,N+3);
    
    Vq2 = Vandermonde1D(N,rq2)/V;
    M = Vq'*diag(wq)*Vq;
    Lq = M\Vf';
    
    Pq = M\(Vq'*diag(wq));
    xq =  Vandermonde1D(N,rq)/Vandermonde1D(N,JacobiGL(0,0,N))*x;
    Pq(abs(Pq)<1e-8) = 0;
    
    VqPq = Vq*Pq;
    
    rp = linspace(-1,1,50); F = eye(N+1);
    Vp = Vandermonde1D(N,rp)/V;
    xp = Vandermonde1D(N,rp)/Vandermonde1D(N,JacobiGL(0,0,N))*x;
    
    global gamma
    gamma = 1.4;
    rhoe = @(rho,m,E) E - .5*m.^2./rho;
    s = @(rho,m,E) log((gamma-1).*rhoe(rho,m,E)./(rho.^gamma));
    
    global V1 V2 V3
    V1 = @(rho,m,E) (-E + rhoe(rho,m,E).*(gamma + 1 - s(rho,m,E)))./(rhoe(rho,m,E));
    V2 = @(rho,m,E) (m)./(rhoe(rho,m,E));
    V3 = @(rho,m,E) (-rho)./(rhoe(rho,m,E));
    
    sV = @(V1,V2,V3) gamma - V1 + V2.^2./(2*V3);
    rhoeV  = @(V1,V2,V3) ((gamma-1)./((-V3).^gamma)).^(1/(gamma-1)).*exp(-sV(V1,V2,V3)/(gamma-1));
    U1 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(-V3);
    U2 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(V2);
    U3 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(1-V2.^2./(2*V3));
    
    k = 1;
    rho = 2 + exp(xq/2).*sin(k*pi*xq);
    rhou = sin(k*pi*xq);
    E = 2 + .5*(rhou.^2)./rho;
    
    rho = VqPq*rho;
    rhou = VqPq*rhou;
    E = VqPq*E;
    
    q1 = Vq2*Pq*V1(rho,rhou,E);
    q2 = Vq2*Pq*V2(rho,rhou,E);
    q3 = Vq2*Pq*V3(rho,rhou,E);
    
    rho = Vq2*Pq*rho;
    rhou = Vq2*Pq*rhou;
    E = Vq2*Pq*E;
    
    rhoV = U1(q1,q2,q3);
    rhouV = U2(q1,q2,q3);
    EV = U3(q1,q2,q3);  
    
    wJq2 = diag(wq2)*(Vq2*J);
    err2 = (rhoV - rho).^2 + (rhouV - rhou).^2 + (EV - E).^2;
    L2err(sk) = sqrt(sum(sum(wJq2.*err2)));
    
    sk = sk + 1;
end
% disp(sprintf('N = %d',N))
print_pgf_coordinates(h,L2err)
C = [h(:).^0 log(h(:))]\log(L2err(:));
C(2)
C2 = diff(log(L2err(end-1:end)))/log(.5)

end
return
loglog(h,L2err,'bo--')
hold on
loglog(h,10*h.^(N+1),'rs-')
C = [h(:).^0 log(h(:))]\log(L2err(:));
C(2)
C2 = diff(log(L2err(end-1:end)))/log(.5)

