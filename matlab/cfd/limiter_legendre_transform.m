clear
N = 6;
r = JacobiGQ(0,0,N);
[rq wq] = JacobiGL(0,0,N+10);
Nq = length(rq);

V = Vandermonde1D(N,r);
Vq = Vandermonde1D(N,rq)/V;
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));

opt = 1;

% entropy variables
if opt==1
    gamma = 1.4;
    rhoe = @(rho,m,E) E - .5*m.^2./rho;
    s = @(rho,m,E) log((gamma-1).*rhoe(rho,m,E)./(rho.^gamma));
    
    V1 = @(rho,m,E) (-E + rhoe(rho,m,E).*(gamma + 1 - s(rho,m,E)))./(rhoe(rho,m,E));
    V2 = @(rho,m,E) (m)./(rhoe(rho,m,E));
    V3 = @(rho,m,E) (-rho)./(rhoe(rho,m,E));
    
    sV = @(V1,V2,V3) gamma - V1 + V2.^2./(2*V3);
    rhoeV  = @(V1,V2,V3) ((gamma-1)./((-V3).^gamma)).^(1/(gamma-1)).*exp(-sV(V1,V2,V3)/(gamma-1));
    U1 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(-V3);
    U2 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(V2);
    U3 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(1-V2.^2./(2*V3));
    
    S = @(rho,m,E) -rho.*s(rho,m,E);
    
    UV = @(V) [U1(V(:,1),V(:,2),V(:,3)) U2(V(:,1),V(:,2),V(:,3)) U3(V(:,1),V(:,2),V(:,3))];
    VU = @(U) [V1(U(:,1),U(:,2),U(:,3)) V2(U(:,1),U(:,2),U(:,3)) V3(U(:,1),U(:,2),U(:,3))];
    SU = @(U) S(U(:,1),U(:,2),U(:,3));
    
    % vector versions
    SV = @(V) SU(UV(V));
    
    % conjugate entropy - Legendre transform
    TV = @(V) sum(UV(V).*V,2) - SV(V);
    TU = @(U) sum(U.*VU(U),2) - SU(U);
    
    rho = 1 + .99*sin(pi*rq);
    u = 1 + .99*sin(pi*rq);
    p = 1 + .99*sin(pi*rq);
    
    m = rho.*u;
    E = p/(gamma-1) + .5*rho.*u.^2;
    
    U = [rho m E];
    
    plot(rq,SU(U),'o--')
    hold on
    plot(rq,TV(VU(U)),'x--')    
    plot(rq,TU(U),'^--')    
    
elseif opt==2
    
    g = 1;
    V1 = @(h,hu) g*h - .5*(hu./h).^2;
    V2 = @(h,hu) hu./h;
    
    U1 = @(v1,v2) (v1+.5*v2.^2)/g;
    U2 = @(v1,v2) U1(v1,v2).*v2;
    
    S = @(h,hu) .5*(hu.^2./h + g*h.^2);
    UV = @(V) [U1(V(:,1),V(:,2)) U2(V(:,1),V(:,2))];
    VU = @(U) [V1(U(:,1),U(:,2)) V2(U(:,1),U(:,2))];
    SU = @(U) S(U(:,1),U(:,2));
    
    % vector versions
    SV = @(V) SU(UV(V));
    
    % conjugate entropy - Legendre transform
    TV = @(V) sum(UV(V).*V,2) - SV(V);
    TU = @(U) sum(U.*VU(U),2) - SU(U);
    
    h = 1 + .5*sin(pi*rq);
    hu = 1 + .5*sin(pi*rq);
    
    U = [h hu];
    
    plot(rq,SU(U),'o--')
    hold on
    plot(rq,TV(VU(U)),'x--')
    plot(rq,.5*h.^2,'^--')
    
end






