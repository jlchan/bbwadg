clear
Globals1D;

CFL = .125/4;
% CFL = .125;
N = 4;
K1D = 40;
FinalTime = 1.8;
% FinalTime = .2;
opt = 3;

global tau
tau = 1;

r = JacobiGL(0,0,N);
% r = JacobiGQ(0,0,N);

[rq wq] = JacobiGL(0,0,N); rq = .9999999999*rq;
[rq wq] = JacobiGQ(0,0,N+1);

% % include boundary nodes for extraction
rq = [-(1-1e-11);rq;(1-1e-11)];
wq = [0;wq;0];
Nq = length(rq);

% evaluate error
rq2 = rq; wq2 = wq;
% [rq2 wq2] = JacobiGQ(0,0,N+4);
% [rq2 wq2] = JacobiGL(0,0,N);

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

global mapP mapB
mapM = reshape(1:Nfp*Nfaces*K,Nfp*Nfaces,K);
mapP = mapM;
for e = 1:K
    mapP(:,e) = [(Nfp*Nfaces)*e-2, Nfp*Nfaces*(e+1)-1];
end
mapP(1) = 1; mapP(end) = Nfp*Nfaces*K;
mapP = reshape(mapP,Nfp*Nfaces,K);
mapB = [1; Nfp*Nfaces*K];

global VqPq tf Dq Vq Pq Dfq
V = Vandermonde1D(N,r);
Dr = GradVandermonde1D(N,r)/V;
Vq = Vandermonde1D(N,rq)/V;
Vf = Vandermonde1D(N,[-1 1])/V;
Vq2 = Vandermonde1D(N,rq2)/V;

M = Vq'*diag(wq)*Vq;
LIFT = M\Vf';

Pq = M\(Vq'*diag(wq));
xq =  Vandermonde1D(N,rq)/Vandermonde1D(N,JacobiGL(0,0,N))*x;
Pq(abs(Pq)<1e-8) = 0;

VqPq = Vq*Pq;
tf = zeros(2,length(rq)); tf(1) = 1; tf(end) = 1;
Dq = Vq*Dr*Pq;
Dfq = .5*Vq*LIFT*diag([-1,1])*(tf - tf*VqPq);

rp = linspace(-1,1,50); F = eye(N+1);
Vp = Vandermonde1D(N,rp)/V;
xp = Vandermonde1D(N,rp)/Vandermonde1D(N,JacobiGL(0,0,N))*x;

rp0 = 0; F0 = V*diag([1;zeros(N,1)])/V;
Vp0 = Vandermonde1D(N,rp0)/V;
xp0 = Vandermonde1D(N,rp0)/Vandermonde1D(N,JacobiGL(0,0,N))*x;

% Vp = Vp0; F = F0; xp = xp0;
% F0 = Vp0*F0;

%% fluxes

global gamma
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

global f1 f2 f3
pavg = @(rhoL,rhoR,uL,uR,EL,ER) avg(rhoL,rhoR)./(2*avg(beta(rhoL,uL,EL),beta(rhoR,uR,ER)));
f1 = @(rhoL,rhoR,uL,uR,EL,ER) logmean(rhoL,rhoR).*avg(uL,uR);
f2 = @(rhoL,rhoR,uL,uR,EL,ER) pavg(rhoL,rhoR,uL,uR,EL,ER) + avg(uL,uR).*f1(rhoL,rhoR,uL,uR,EL,ER);
f3 = @(rhoL,rhoR,uL,uR,EL,ER) f1(rhoL,rhoR,uL,uR,EL,ER)...
    .*(1./(2*(gamma-1).*logmean(beta(rhoL,uL,EL),beta(rhoR,uR,ER))) - .5*avg(uL.^2,uR.^2)) ...
    + avg(uL,uR).*f2(rhoL,rhoR,uL,uR,EL,ER);


% Jacobian on time derivative
dVdU11 = @(rho,m,E)(1.0./(E.*rho.*2.0-m.^2).^2.*(gamma.*m.^4+m.^4+E.^2.*gamma.*rho.^2.*4.0-E.*gamma.*m.^2.*rho.*4.0))./rho;
dVdU12 = @(rho,m,E)m.^3.*1.0./(E.*rho.*2.0-m.^2).^2.*-2.0;
dVdU13 = @(rho,m,E)rho.*(E.*rho-m.^2).*1.0./(E.*rho.*2.0-m.^2).^2.*-4.0;
dVdU22 = @(rho,m,E)rho.*(E.*rho.*2.0+m.^2).*1.0./(E.*rho.*2.0-m.^2).^2.*2.0;
dVdU23 = @(rho,m,E)m.*rho.^2.*1.0./(E.*rho.*2.0-m.^2).^2.*-4.0;
dVdU33 = @(rho,m,E)rho.^3.*1.0./(E.*rho.*2.0-m.^2).^2.*4.0;

% dUdV11 = @(V1,V2,V3)-(V3.*exp(-(-V1+gamma+(V2.^2.*(1.0./2.0))./V3)./(gamma-1.0)).*((-V3).^(-gamma).*(gamma-1.0)).^(1.0./(gamma-1.0)))./(gamma-1.0);
% dUdV12 = @(V1,V2,V3)(V2.*exp(-(-V1+gamma+(V2.^2.*(1.0./2.0))./V3)./(gamma-1.0)).*((-V3).^(-gamma).*(gamma-1.0)).^(1.0./(gamma-1.0)))./(gamma-1.0);
% dUdV13 = @(V1,V2,V3)(exp((V1.*V3-V3.*gamma-V2.^2.*(1.0./2.0))./(V3.*(gamma-1.0))).*((-V3).^(-gamma).*(gamma-1.0)).^(1.0./(gamma-1.0)).*(V3.*2.0-V2.^2).*(1.0./2.0))./(V3.*(gamma-1.0));
% dUdV22 = @(V1,V2,V3)-(exp((V1.*V3-V3.*gamma-V2.^2.*(1.0./2.0))./(V3.*(gamma-1.0))).*((-V3).^(-gamma).*(gamma-1.0)).^(1.0./(gamma-1.0)).*(V3-V3.*gamma+V2.^2))./(V3.*(gamma-1.0));
% dUdV23 = @(V1,V2,V3)(V2.*1.0./V3.^2.*exp((V1.*V3-V3.*gamma-V2.^2.*(1.0./2.0))./(V3.*(gamma-1.0))).*((-V3).^(-gamma).*(gamma-1.0)).^(1.0./(gamma-1.0)).*(V3.*gamma.*2.0-V2.^2).*(-1.0./2.0))./(gamma-1.0);
% dUdV33 = @(V1,V2,V3)(1.0./V3.^3.*exp((V1.*V3-V3.*gamma-V2.^2.*(1.0./2.0))./(V3.*(gamma-1.0))).*((-V3).^(-gamma).*(gamma-1.0)).^(1.0./(gamma-1.0)).*(V3.^2.*gamma.*4.0+V2.^4-V2.^2.*V3.*gamma.*4.0).*(-1.0./4.0))./(gamma-1.0);


%%
res1 = 0;
res2 = 0;
res3 = 0;

CN = (N+1)^2/2;
h = 2/K1D;
dt = CFL * h/CN;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

global rhoL rhoR mL mR EL ER
global useBC; useBC = 1;

if opt==0
    
    % sine solution
    t = 0;
    rhoex = @(x) (2 + sin(pi*(x - t)));
    uex = @(x) ones(size(x));
    pex = @(x) ones(size(x));
    mex = @(x) rhoex(x).*uex(x);
    
    Eex = @(x) pex(x)./(gamma-1) + .5*rhoex(x).*uex(x).^2;
    
    useBC = 0;
    mapP(1) = mapM(end); mapP(end) = mapM(1);
    vmapP(1) = vmapM(end); vmapP(end) = vmapM(1);

elseif opt==1
    
    % pulse condition
    x0 = 0;
    rhoex = @(x) (2 + (abs(x-x0) < .5));
    pex = @(x) rhoex(x).^gamma;
    uex = @(x) zeros(size(x));
    Eex = @(x) pex(x)/(gamma-1) + .5*rhoex(x).*uex(x).^2;
    mex = @(x) rhoex(x).*uex(x);
    
    useBC = 0;
    mapP(1) = mapM(end); mapP(end) = mapM(1);
    vmapP(1) = vmapM(end); vmapP(end) = vmapM(1);
    
elseif opt==2
    
    % sod initial condition
    rhoex = @(x) (1*(x < 0) + .125*(x > 0));
    pex = @(x) (1*(x < 0) + .1*(x > 0));
    uex = @(x) zeros(size(x));
    Eex = @(x) pex(x)/(gamma-1) + .5*rhoex(x).*uex(x).^2;
    mex = @(x) rhoex(x).*uex(x);
    
    % BCs
    rhoL = 1; rhoR = .125;
    pL = 1; pR = .1;
    uL = 0; uR = 0;
    mL = 0; mR = 0;
    EL = pL/(gamma-1) + .5*rhoL*uL.^2;
    ER = pR/(gamma-1) + .5*rhoR*uR.^2;

elseif opt==3
    
    rhoRex = @(x) 1 + .2*sin(5*x);
    rhoL = 3.857143; rhoR = rhoRex(x(end));
    uL   = 2.629369; uR = 0;
    pL   = 10.3333;  pR = 1;

% rhoL = 4; rhoR = rhoRex(x(end));
% uL   = 2; uR = 0;
% pL   = 8;  pR = 1;
        
    rhoex = @(x) rhoL*(x < -4) + (rhoRex(x)).*(x >= -4);
    uex = @(x) uL*(x < -4);
    pex = @(x) pL*(x < -4) + pR*(x >= -4);
    
    Eex = @(x) (pex(x)/(gamma-1) + .5*rhoex(x).*uex(x).^2);
    mex = @(x) (rhoex(x).*uex(x));
    
    mL = rhoL*uL; mR = rhoR*uR;
    EL = pL/(gamma-1) + .5*rhoL*uL.^2;
    ER = pR/(gamma-1) + .5*rhoR*uR.^2;
    
end

q1 = Pq*V1(rhoex(xq),mex(xq),Eex(xq));
q2 = Pq*V2(rhoex(xq),mex(xq),Eex(xq));
q3 = Pq*V3(rhoex(xq),mex(xq),Eex(xq));


% plot(xp,U1(Vp*q1,Vp*q2,Vp*q3))
% hold on
% plot(xp,U2(Vp*q1,Vp*q2,Vp*q3))
% plot(xp,U3(Vp*q1,Vp*q2,Vp*q3))
% return
% useBC = 0;
% mapP(1) = mapM(end); mapP(end) = mapM(1);
% vmapP(1) = vmapM(end); vmapP(end) = vmapM(1);


wJq = diag(wq)*(Vq*J);

figure(1)
for i = 1:Nsteps  
    
    if (i < 150)
        CFL = .125/4;
        dt = CFL*h/CN;
        Nsteps = ceil(FinalTime/dt);
        dt = FinalTime/Nsteps;
    elseif i == 150
        CFL = .125/2.5;
        dt = CFL*h/N;        
    end
    if (dt*i > FinalTime)
        dt = FinalTime-(i-1)*dt;        
        i=Nsteps;
    end
    
    for INTRK = 1:5
        
        % interpolate to quadrature
        q1q = Vq*q1;
        q2q = Vq*q2;
        q3q = Vq*q3;
        
        rhoq = U1(q1q,q2q,q3q);
        mq   = U2(q1q,q2q,q3q);
        Eq   = U3(q1q,q2q,q3q);
        
        [rhs1 rhs2 rhs3] = rhsEuler(rhoq,mq,Eq,q1q,q2q,q3q);
        
        % apply WADG
        r1 = Vq*rhs1;
        r2 = Vq*rhs2;
        r3 = Vq*rhs3;
              
%         F = Vq*Pq;
        F = eye(Nq); 
        A11 = F*dVdU11(rhoq,mq,Eq);
        A12 = F*dVdU12(rhoq,mq,Eq);
        A13 = F*dVdU13(rhoq,mq,Eq);
        A22 = F*dVdU22(rhoq,mq,Eq);
        A23 = F*dVdU23(rhoq,mq,Eq);
        A33 = F*dVdU33(rhoq,mq,Eq);
        
        rhs1 = Pq*(A11.*r1 + A12.*r2 + A13.*r3);
        rhs2 = Pq*(A12.*r1 + A22.*r2 + A23.*r3);
        rhs3 = Pq*(A13.*r1 + A23.*r2 + A33.*r3);                        
        
%         keyboard
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        q1 = q1 + rk4b(INTRK)*res1;
        q2 = q2 + rk4b(INTRK)*res2;
        q3 = q3 + rk4b(INTRK)*res3;
        
    end
    
    q1q = Vq*q1;
    q2q = Vq*q2;
    q3q = Vq*q3;    
    rhoq = U1(q1q,q2q,q3q);
    mq   = U2(q1q,q2q,q3q);
    Eq   = U3(q1q,q2q,q3q);    
    sq = s(rhoq,mq,Eq);
    Sq = -rhoq.*sq/(gamma-1);
    energy(i) = sum(sum(wJq.*Sq));
    
    if mod(i,10)==0 || i==Nsteps
        
        q1p = Vp*q1;
        q2p = Vp*q2;
        q3p = Vp*q3;
        
        rhop = U1(q1p,q2p,q3p);
        mp = U2(q1p,q2p,q3p);
        Ep = U3(q1p,q2p,q3p);
        up = mp./rhop;
        pp = (gamma-1)*(Ep-.5*rhop.*up.^2);
        plot(xp,rhop,'b-','linewidth',2)
        hold on
        plot(xp,up,'r-.','linewidth',2)
        plot(xp,pp,'k--','linewidth',2)
        
        q1q = Vq*q1;
        q2q = Vq*q2;
        q3q = Vq*q3;
        rho0 = sum(wJq.*U1(q1q,q2q,q3q),1)./sum(wJq,1);        
        m0 = sum(wJq.*U2(q1q,q2q,q3q),1)./sum(wJq,1);        
        E0 = sum(wJq.*U3(q1q,q2q,q3q),1)./sum(wJq,1);
        u0 = m0./rho0;
        p0 = (gamma-1)*(E0-.5*rho0.*u0.^2);
        plot(xp0,rho0,'bo','linewidth',2)
        plot(xp0,u0,'ro','linewidth',2)
        plot(xp0,p0,'ko','linewidth',2)
        
        title(sprintf('Time = %f, step = %d out of %d',dt*i,i,Nsteps))
        %         axis([-5 5 -1 7])
        hold off
        drawnow
    end
end

% sine solution
t = FinalTime;
rhoex = @(x) (2 + sin(pi*(x - t)));
uex = @(x) ones(size(x));
pex = @(x) ones(size(x));
mex = @(x) rhoex(x).*uex(x);    
Eex = @(x) pex(x)./(gamma-1) + .5*rhoex(x).*uex(x).^2;

q1q = Vq*q1;
q2q = Vq*q2;
q3q = Vq*q3;
e1 = rhoex(xq) - U1(q1q,q2q,q3q);
e2 = mex(xq) - U2(q1q,q2q,q3q);
e3 = Eex(xq) - U3(q1q,q2q,q3q);
L2err = sqrt(sum(sum(wJq.*(e1.^2+e2.^2+e3.^2))));
L2err
    

% figure(2)
% semilogy(dt*(1:Nsteps),abs(energy),'--')
% hold on


return

function [rhs1 rhs2 rhs3] = rhsEuler(rhoq,mq,Eq,q1q,q2q,q3q)

Globals1D
global rhoL rhoR mL mR EL ER gamma
global mapP mapB
global VqPq tf Dq Vq Pq Dfq
global f1 f2 f3

uq = mq./rhoq;
pq = (gamma-1)*(Eq-.5*rhoq.*uq.^2);

Nq = size(rhoq,1);
rhoM = rhoq([1 Nq],:); uM = uq([1 Nq],:); EM = Eq([1 Nq],:);
mM = rhoM.*uM;
rhoP = reshape(rhoM(mapP),Nfp*Nfaces,K);
uP = reshape(uM(mapP),Nfp*Nfaces,K);
EP = reshape(EM(mapP),Nfp*Nfaces,K);
mP = rhoP.*uP;
pM = pq([1 Nq],:);


q1M = q1q([1 Nq],:); 
q2M = q2q([1 Nq],:); 
q3M = q3q([1 Nq],:);
q1P = reshape(q1M(mapP),Nfp*Nfaces,K);
q2P = reshape(q2M(mapP),Nfp*Nfaces,K);
q3P = reshape(q3M(mapP),Nfp*Nfaces,K);


global useBC
if useBC
    rhoP(mapB) = [rhoL rhoR];
    uP(mapB) = [mL/rhoL mR/rhoR];
    mP(mapB) = [mL mR];
    EP(mapB) = [EL ER];
end

% compute fluxes
f1f = f1(rhoM,rhoP,uM,uP,EM,EP);
f2f = f2(rhoM,rhoP,uM,uP,EM,EP);
f3f = f3(rhoM,rhoP,uM,uP,EM,EP);

rhs1 = zeros(Np,K);
rhs2 = zeros(Np,K);
rhs3 = zeros(Np,K);
for e = 1:K
    [rhox rhoy] = meshgrid(rhoq(:,e));
    [ux uy] = meshgrid(uq(:,e));
    [Ex Ey] = meshgrid(Eq(:,e));
    
    [rhoxf rhoyf] = meshgrid(rhoq(:,e),tf*rhoq(:,e));
    [uxf uyf] = meshgrid(uq(:,e),tf*uq(:,e));
    [Exf Eyf] = meshgrid(Eq(:,e),tf*Eq(:,e));
    
    flux1 = .5*diag([-1 1])*f1f(:,e);
    flux2 = sum((.5*diag([-1 1])*(-tf*VqPq)).*f1(rhoxf,rhoyf,uxf,uyf,Exf,Eyf),2);
    rhs1(:,e) = rx(:,e).*(Pq*sum(Dq.*f1(rhox,rhoy,ux,uy,Ex,Ey),2)) + Fscale(1,e)*(Pq*sum(Dfq.*f1(rhox,rhoy,ux,uy,Ex,Ey),2)) + LIFT*(Fscale(:,e).*(flux1+flux2));
    
    flux1 = .5*diag([-1 1])*f2f(:,e);
    flux2 = sum((.5*diag([-1 1])*(-tf*VqPq)).*f2(rhoxf,rhoyf,uxf,uyf,Exf,Eyf),2);
    rhs2(:,e) = rx(:,e).*(Pq*sum(Dq.*f2(rhox,rhoy,ux,uy,Ex,Ey),2)) + Fscale(1,e)*(Pq*sum(Dfq.*f2(rhox,rhoy,ux,uy,Ex,Ey),2)) + LIFT*(Fscale(:,e).*(flux1+flux2));
    
    flux1 = .5*diag([-1 1])*f3f(:,e);
    flux2 = sum((.5*diag([-1 1])*(-tf*VqPq)).*f3(rhoxf,rhoyf,uxf,uyf,Exf,Eyf),2);
    rhs3(:,e) = rx(:,e).*(Pq*sum(Dq.*f3(rhox,rhoy,ux,uy,Ex,Ey),2)) + Fscale(1,e)*(Pq*sum(Dfq.*f3(rhox,rhoy,ux,uy,Ex,Ey),2)) + LIFT*(Fscale(:,e).*(flux1+flux2));
    
end

global tau

% local lax penalty
cvel = sqrt(gamma*pM./rhoM);
lm   = (abs(uM) + cvel);
Lfc  = max(lm(mapP),lm);
% Lfc  = .5*(lm(mapP)+lm);

rhs1 = -2*rhs1 + .5*tau*LIFT*(Fscale.*Lfc.*(rhoP-rhoM));
rhs2 = -2*rhs2 + .5*tau*LIFT*(Fscale.*Lfc.*(mP-mM));
rhs3 = -2*rhs3 + .5*tau*LIFT*(Fscale.*Lfc.*(EP-EM));
% rhs1 = -2*rhs1 + .5*tau*LIFT*(Fscale.*Lfc.*(q1P-q1M));
% rhs2 = -2*rhs2 + .5*tau*LIFT*(Fscale.*Lfc.*(q2P-q2M));
% rhs3 = -2*rhs3 + .5*tau*LIFT*(Fscale.*Lfc.*(q3P-q3M));

end