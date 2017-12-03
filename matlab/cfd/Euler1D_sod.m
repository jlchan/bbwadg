clear
Globals1D;

projectV = 1;
CFL = .125;
N = 4;
K1D = 16;
tau = 1/2;
FinalTime = .2;

r = JacobiGL(0,0,N);
% r = JacobiGQ(0,0,N);

[rq wq] = JacobiGL(0,0,N);
[rq wq] = JacobiGQ(0,0,N+1);
% rq = .999*rq;

% % include boundary nodes for extraction
rq = [-1;rq;1];
wq = [0;wq;0];
Nq = length(rq);

% evaluate error
rq2 = rq; wq2 = wq;
% [rq2 wq2] = JacobiGQ(0,0,N+4);
% [rq2 wq2] = JacobiGL(0,0,N);

[Nv, VX, K, EToV] = MeshGen1D(-.5,.5,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

mapM = reshape(1:Nfp*Nfaces*K,Nfp*Nfaces,K);
mapP = mapM;
for e = 1:K
    mapP(:,e) = [(Nfp*Nfaces)*e-2, Nfp*Nfaces*(e+1)-1];
end
mapP(1) = 1; mapP(end) = Nfp*Nfaces*K;
% mapP(1) = mapM(end); mapP(end) = mapM(1);
% vmapP(1) = vmapM(end); vmapP(end) = vmapM(1);
mapB = [1; Nfp*Nfaces*K];

V = Vandermonde1D(N,r);
Dr = GradVandermonde1D(N,r)/V;
Vq = Vandermonde1D(N,rq)/V;
Vf = Vandermonde1D(N,[-1 1])/V;
Vq2 = Vandermonde1D(N,rq2)/V;

M = Vq'*diag(wq)*Vq;
LIFT = M\Vf';

Vrq = GradVandermonde1D(N,rq)/V;
Pq = M\(Vq'*diag(wq));
Prq = M\(Vrq'*diag(wq));
xq =  Vandermonde1D(N,rq)/Vandermonde1D(N,JacobiGL(0,0,N))*x;
Pq(abs(Pq)<1e-8) = 0;

nxJ = nx.*Fscale;
du = @(uf) reshape(uf(mapP) - uf(mapM), Nfp*Nfaces, K);
Dh = @(u) rx.*(Dr*u) + .5*LIFT*(du(Vf*u).*nxJ);

D = zeros((N+1)*K);
uu = zeros(N+1,K);
for i = 1:(N+1)*K
    uu(i) = 1;
    tmp = Dh(uu);
    D(:,i) = tmp(:);
    uu(i) = 0;
end
D(abs(D)<1e-8) = 0;

rp = linspace(-1,1,100)';
Vp = Vandermonde1D(N,rp)/V;
xp = Vandermonde1D(N,rp)/Vandermonde1D(N,JacobiGL(0,0,N))*x;

Dx = kron(spdiag(rx(1,:)),Dr);
Dq = kron(speye(K),Vq)*D*kron(speye(K),Pq);
Dq(abs(Dq)<1e-8) = 0;
Dq = sparse(Dq);
PqG = kron(speye(K),Pq);
MM = kron(spdiag(J(1,:)),Vq'*diag(wq)*Vq);

VfG = kron(speye(K),Vf);
VqG = kron(speye(K),Vq);
VpG = kron(speye(K),Vp);
xq = xq(:); xp = xp(:);
wJq = diag(wq)*(Vq*J); wJq = wJq(:);
wJq2 = diag(wq2)*(Vq2*J); wJq2 = wJq2(:);
Vq2G = kron(speye(K),Vq2);
xq2 = Vq2G*x(:);

W = diag(wJq);

% chen/shu ops
% global lift
L = zeros((N+1)*K,Nfp*Nfaces*K);
uu = zeros(Nfp*Nfaces,K);
for i = 1:Nfp*Nfaces*K
    uu(i) = 1;
    tmp = LIFT*(reshape(uu,Nfp*Nfaces,K).*nx.*Fscale);
    L(:,i) = tmp(:);
    uu(i) = 0;
end
L(abs(L)<1e-8) = 0; L = sparse(L);

S = zeros((N+1)*K,Nfp*Nfaces*K);
uu = zeros(Nfp*Nfaces,K);
for i = 1:Nfp*Nfaces*K
    uu(i) = 1;
    tmp = LIFT*(reshape(uu,Nfp*Nfaces,K).*Fscale);
    S(:,i) = tmp(:);
    uu(i) = 0;
end
S(abs(S)<1e-8) = 0; S = sparse(S);

% extract from all points
tM = zeros(Nfp*Nfaces*K,Nq*K);
tP = zeros(Nfp*Nfaces*K,Nq*K);
uu = zeros(Nq,K);
for i = 1:Nq*K
    uu(i) = 1;
    uf = uu([1; Nq],:);
    
    tmp = uf(mapM);
    tM(:,i) = tmp(:);
    tmp = uf(mapP);
    tP(:,i) = tmp(:);
    
    uu(i) = 0;
end

nx = nx(:);
%Dh = .5*(VqG*Dx*PqG - VqG*(MM\(Dx'*VqG'*W))) + .5*VqG*L*(tM);
Dh = VqG*Dx*PqG + .5*VqG*L*(tM - tM*VqG*PqG);
Fh = .5*tM'*diag(nx(:))*(tP - tM*VqG*PqG);
Fh2 = .5*tM'*diag(nx(:))*(-tM*VqG*PqG);
Pfh = MM\(VqG');
% Fh = .5*diag(nx(:))*(tP - tM*VqG*PqG);
% Pfh = MM\(VfG');

% Dh = VqG*Dx*PqG + .5*VqG*L*(tP-tM)*VqG*PqG;
% Fh = 0*Fh;

%% fluxes

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

avg = @(uL,uR) .5*(uL+uR);
pfun = @(rho,u,E) (gamma-1)*(E-.5*rho.*u.^2);
beta = @(rho,u,E) rho./(2*pfun(rho,u,E));
f1 = @(rhoL,rhoR,uL,uR,EL,ER) logmean(rhoL,rhoR).*avg(uL,uR);
f2 = @(rhoL,rhoR,uL,uR,EL,ER) avg(rhoL,rhoR)./(2*avg(beta(rhoL,uL,EL),beta(rhoR,uR,ER))) + avg(uL,uR).*f1(rhoL,rhoR,uL,uR,EL,ER);
f3 = @(rhoL,rhoR,uL,uR,EL,ER) f1(rhoL,rhoR,uL,uR,EL,ER)...
    .*(1./(2*(gamma-1).*logmean(beta(rhoL,uL,EL),beta(rhoR,uR,ER))) - .5*avg(uL.^2,uR.^2)) ...
    + avg(uL,uR).*f2(rhoL,rhoR,uL,uR,EL,ER);


%%
res1 = 0;
res2 = 0;
res3 = 0;

CN = (N+1)^2/2;
h = 2/K1D;
dt = CFL * h/CN;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

opt = 1;
if opt==1 % sod
    
    rho = PqG*(1*(xq < 0) + .125*(xq > 0));
    p = PqG*(1*(xq < 0) + .1*(xq > 0));
    u = 0*rho;
    E = p/(gamma-1) + .5*rho.*u.^2;
    m = rho.*u;
    
    rhoL = 1; rhoR = .125;
    pL = 1; pR = .1;
    uL = 0; uR = 0;
    mL = 0; mR = 0;
    EL = pL/(gamma-1) + .5*rhoL*uL.^2;
    ER = pR/(gamma-1) + .5*rhoR*uR.^2;

end


F = kron(eye(K),V*diag([ones(N-1,1);1;1])/V);
F0 = kron(eye(K),V*diag([1;zeros(N,1)])/V);

% x = x(:);
% plot(x,rho)
% hold on
% plot(x,m)
% plot(x,E)
% return

energy = zeros(Nsteps,1);
figure(1)
for i = 1:Nsteps    
    for INTRK = 1:5
        
        % interpolate to quadrature
        rhoq = VqG*rho;
        mq = VqG*m;
        Eq = VqG*E;
        
        % project entropy variables
        q1 = VqG*PqG*V1(rhoq,mq,Eq);
        q2 = VqG*PqG*V2(rhoq,mq,Eq);
        q3 = VqG*PqG*V3(rhoq,mq,Eq);
        if projectV
            rhoq = U1(q1,q2,q3);
            mq   = U2(q1,q2,q3);
            Eq   = U3(q1,q2,q3);
        end
        uq = mq./rhoq;
        pq = (gamma-1)*(Eq-.5*rhoq.*uq.^2);
        
        [rhox rhoy] = meshgrid(rhoq);
        [ux uy] = meshgrid(uq);
        [Ex Ey] = meshgrid(Eq);
        
        % local lax penalty
        cvel = sqrt(gamma*pq./rhoq);
        lm   = (abs(uq) + cvel);
        Lfc   = max(tM*lm,tP*lm);
        
        rhoM = tM*rhoq; uM = tM*uq; EM = tM*Eq;                
        rhoP = tP*rhoq; uP = tP*uq; mP = tP*mq; 
        EP = tP*Eq; 
        rhoP(mapB) = [rhoL rhoR]; 
        uP(mapB) = [uL uR]; 
        mP(mapB) = [mL mR]; 
        EP(mapB) = [EL ER];        
        
        flux1 = .5*tM'*diag(nx(:))*(f1(rhoM,rhoP,uM,uP,EM,EP));
        flux2 = sum(Fh2.*f1(rhox,rhoy,ux,uy,Ex,Ey),2);
        d1 = PqG*(sum(Dh.*f1(rhox,rhoy,ux,uy,Ex,Ey),2)) + Pfh*(flux1 + flux2);
        rhs1 = -2*d1 + tau*S*(Lfc.*(rhoP-tM*rhoq));                       
        
        flux1 = .5*tM'*diag(nx(:))*(f2(rhoM,rhoP,uM,uP,EM,EP));
        flux2 = sum(Fh2.*f2(rhox,rhoy,ux,uy,Ex,Ey),2);
        d2 = PqG*(sum(Dh.*f2(rhox,rhoy,ux,uy,Ex,Ey),2)) + Pfh*(flux1 + flux2);        
        rhs2 = -2*d2 + tau*S*(Lfc.*(mP-tM*mq));
        
        flux1 = .5*tM'*diag(nx(:))*(f3(rhoM,rhoP,uM,uP,EM,EP));
        flux2 = sum(Fh2.*f3(rhox,rhoy,ux,uy,Ex,Ey),2);        
        d3 = PqG*(sum(Dh.*f3(rhox,rhoy,ux,uy,Ex,Ey),2)) + Pfh*(flux1 + flux2);
        rhs3 = -2*d3 + tau*S*(Lfc.*(EP-tM*Eq));                

        if (INTRK==5)
            rhstest(i) = sum(wJq.*(q1.*(VqG*rhs1) + q2.*(VqG*rhs2) + q3.*(VqG*rhs3)));
        end
        
%         if INTRK==0
%             Dh = 0*VqG*Dx*PqG + .5*VqG*L*(tM - tM*VqG*PqG);
%             [PqG*(sum(Dh.*f1(rhox,rhoy,ux,uy,Ex,Ey),2)) PqG*(sum(Dh.*f2(rhox,rhoy,ux,uy,Ex,Ey),2)) PqG*(sum(Dh.*f3(rhox,rhoy,ux,uy,Ex,Ey),2))]
% %             [rhs1 rhs2 rhs3]
% %             keyboard
%         end
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        rho = rho + rk4b(INTRK)*res1;
        m = m + rk4b(INTRK)*res2;
        E = E + rk4b(INTRK)*res3;
        
    end;
    
    rhoq = Vq2G*rho;
    mq = Vq2G*m;
    Eq = Vq2G*E;
    uq = mq./rhoq;
    sq = s(rhoq,mq,Eq);
    %     sq = log(pfun(rhoq,uq,Eq)) - gamma*log(rhoq);
    Sq = -rhoq.*sq/(gamma-1);
    energy(i) = (sum(wJq2.*Sq));
    
    if mod(i,10)==0 || i==Nsteps
        rhop = VpG*F*rho;
        mp = VpG*F*m;
        plot(xp,rhop,'o')
        hold on
        plot(xp,mp./rhop,'x')
        title(sprintf('Time = %f',dt*i))
%         axis([-1/2 1/2 -1 4])
axis([-.5 .5 0 1.5])
        
        hold off
        drawnow
        %         pause
    end
end

max(abs(rhstest))
