clear
clear -global
Globals1D;

projectV = 1;
CFL = .75;
% CFL = .125/2.5;
N = 1;
K1D = 64;

FinalTime = 1.8;
FinalTime = .2;
FinalTime = 5;
opt = 3;

global tau
tau = 1;

% r = JacobiGL(0,0,N);
r = JacobiGL(0,0,N);

[rq wq] = JacobiGL(0,0,N); rq = rq*(1-1e-11);
% [rq wq] = JacobiGL(0,0,N+1); rq = rq*(1-1e-11);
[rq wq] = JacobiGQ(0,0,N+1);
[rqref wqref] = JacobiGQ(0,0,N+1);

% % check negative weights
% rq = linspace(-1,1,N+4)';
% wq = wqref'*(Vandermonde1D(length(rq)-1,rqref)/Vandermonde1D(length(rq)-1,rq)); wq = wq(:);

Nq = length(rq);

if opt==1 || opt==0
    [Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
elseif opt==2
    [Nv, VX, K, EToV] = MeshGen1D(-.5,.5,K1D);
elseif opt==3
    [Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
elseif opt==4
    [Nv, VX, K, EToV] = MeshGen1D(0,1,K1D);
end

% VX = VX*5;

% Initialize solver and construct grid and metric
StartUp1D;

global VqPq Ef Vq Pq Lq
global Vf wq 
V = Vandermonde1D(N,r);
Dr = GradVandermonde1D(N,r)/V;
Vq = Vandermonde1D(N,rq)/V;
Vf = Vandermonde1D(N,[-1 1])/V;

M = Vq'*diag(wq)*Vq;
Lq = M\Vf';

Pq = M\(Vq'*diag(wq));
xq =  Vandermonde1D(N,rq)/Vandermonde1D(N,JacobiGL(0,0,N))*x;
Pq(abs(Pq)<1e-8) = 0;

VqPq = Vq*Pq;
VqLq = Vq*Lq;
VfPq = Vf*Pq;
Drq = Vq*Dr*Pq;

global DNr
nrJ = [-1;1];

Qr = diag(wq)*Vq*Dr*Pq;
Qrskew = (Qr-Qr');
QNr = .5*[Qrskew (Vf*Pq)'*diag([-1;1])
    -diag([-1;1])*Vf*Pq zeros(2)];
WN = diag([wq;1;1]);
PN = M\[Vq' Vf']*WN;
DNr = WN\QNr;

W = diag([wq;1;1]);

Ef = zeros(2,length(rq)); Ef(1) = 1; Ef(end) = 1;

global Vp xp
rp = linspace(-1,1,25); 
Vp = Vandermonde1D(N,rp)/V;
xp = Vandermonde1D(N,rp)/Vandermonde1D(N,JacobiGL(0,0,N))*x;

% filtering
global F
d = ones(N+1,1);
% d(end) = 0;
F = V*(diag(d)/V);


%% maps

global mapP mapB
mapM = reshape(1:Nfp*Nfaces*K,Nfp*Nfaces,K);
mapP = mapM;
for e = 1:K
    mapP(:,e) = [(Nfp*Nfaces)*e-2, Nfp*Nfaces*(e+1)-1];
end
mapP(1) = 1; mapP(end) = Nfp*Nfaces*K;
mapP = reshape(mapP,Nfp*Nfaces,K);
mapB = [1; Nfp*Nfaces*K];


%% fluxes

global gamma
gamma = 1.4;
rhoe = @(rho,m,E) E - .5*m.^2./rho;
s = @(rho,m,E) log((gamma-1).*rhoe(rho,m,E)./(rho.^gamma));

global U1 U2 U3 V1 V2 V3
V1 = @(rho,m,E) (-E + rhoe(rho,m,E).*(gamma + 1 - s(rho,m,E)))./(rhoe(rho,m,E));
V2 = @(rho,m,E) (m)./(rhoe(rho,m,E));
V3 = @(rho,m,E) (-rho)./(rhoe(rho,m,E));

sV = @(V1,V2,V3) gamma - V1 + V2.^2./(2*V3);
rhoeV  = @(V1,V2,V3) ((gamma-1)./((-V3).^gamma)).^(1/(gamma-1)).*exp(-sV(V1,V2,V3)/(gamma-1));
U1 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(-V3);
U2 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(V2);
U3 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(1-V2.^2./(2*V3));

global avg pfun beta
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


%% set initial cond

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
    
    rhoq = rhoex(xq);
    uq = uex(xq);
    pq = pex(xq);
    mq = mex(xq);
    Eq = Eex(xq);
    
    useBC = 0;
    mapP(1) = mapM(end); mapP(end) = mapM(1);
    vmapP(1) = vmapM(end); vmapP(end) = vmapM(1);
    
elseif opt==1
    
    % pulse condition
    x0 = 0;
    rhoq = (2 + (abs(xq-x0) < .5));
    pq = rhoq.^gamma;
    uq = 0*rhoq;
    
    useBC = 0;
    mapP(1) = mapM(end); mapP(end) = mapM(1);
    vmapP(1) = vmapM(end); vmapP(end) = vmapM(1);
    
elseif opt==2
    
    % BCs
    rhoL = 1; rhoR = .125;
    pL = 1; pR = .1;
    uL = 0; uR = 0;
    mL = 0; mR = 0;
    EL = pL/(gamma-1) + .5*rhoL*uL.^2;
    ER = pR/(gamma-1) + .5*rhoR*uR.^2;
    
    % sod initial condition
    rhoq = (rhoL*(xq < 0) + rhoR*(xq > 0));
    pq = (pL*(xq < 0) + pR*(xq > 0));
    uq = 0*xq;
    Eq = pq/(gamma-1) + .5*rhoq.*uq.^2;
    mq = rhoq.*uq;
    
    
elseif opt==3
    useBC = 0;
    mapP(1) = mapM(end); mapP(end) = mapM(1);
    vmapP(1) = vmapM(end); vmapP(end) = vmapM(1);
        
    du = cos(pi*xq);
    rhoq = 2 + sin(pi*xq);
    pq = ones(size(xq));
    uq = .1*(ones(size(xq)) + .5*du);
end

rho = Pq*(rhoq);
u = Pq*(uq);
p = Pq*(pq);
E = Pq*(pq/(gamma-1) + .5*rhoq.*uq.^2);
m = Pq*(rhoq.*uq);

wJq = diag(wq)*((Vandermonde1D(N,rq)/V)*J);

%% entropy functions

global Sfun SU Sref
VV = @(U) [V1(U(:,1),U(:,2),U(:,3)),V2(U(:,1),U(:,2),U(:,3)),V3(U(:,1),U(:,2),U(:,3))];
UU = @(V) [U1(V(:,1),V(:,2),V(:,3)),U2(V(:,1),V(:,2),V(:,3)),U3(V(:,1),V(:,2),V(:,3))];

Sref = max(max(-rhoq.*s(Vq*rho,Vq*m,Vq*E)));
Sfun = @(rho,m,E) -rho.*s(rho,m,E) - Sref - 1;
SU = @(U) Sfun(U(:,1),U(:,2),U(:,3));

%%

% compute reference entropy
rhoq = Vq*rho;
mq = Vq*m;
Eq = Vq*E;

S0 = Sfun([Vq;Vf]*rho,[Vq;Vf]*m,[Vq;Vf]*E);
Sref = max(abs(S0(:)));

res1 = 0;
res2 = 0;
res3 = 0;

CN = (N+1)^2/2;
h = (max(VX)-min(VX))/K1D;
dt = CFL * h/CN;
Nsteps = ceil(FinalTime/dt);


VNPq = [Vq;Vf]*Pq;

interval = 5;
Nsteps = interval*ceil(Nsteps/interval);
dt = FinalTime/Nsteps;
Usnap = zeros(3*(N+1)*K,Nsteps/interval+2);
Vsnap = zeros(3*(N+1)*K,Nsteps/interval+2);

rhoq = Vq*rho;
mq = Vq*m;
Eq = Vq*E;
q1 = Pq*V1(rhoq,mq,Eq);
q2 = Pq*V2(rhoq,mq,Eq);
q3 = Pq*V3(rhoq,mq,Eq);
Usnap(:,1) = [rho(:);m(:);E(:)];
Vsnap(:,1) = [q1(:);q2(:);q3(:)];

figure(1)
sk = 2;
for i = 1:Nsteps        
    
    for INTRK = 1:5
        
        % interpolate to quadrature
        rhoq = Vq*rho;
        mq = Vq*m;
        Eq = Vq*E;
        
        % project entropy variables
        q1 = VNPq*V1(rhoq,mq,Eq);
        q2 = VNPq*V2(rhoq,mq,Eq);
        q3 = VNPq*V3(rhoq,mq,Eq);
        
        rhoq = U1(q1,q2,q3);
        mq   = U2(q1,q2,q3);
        Eq   = U3(q1,q2,q3);
               
        [rhs1 rhs2 rhs3] = rhsEuler(rhoq,mq,Eq,rho,m,E,i,INTRK);
        
        if (INTRK==5)
            rhstest(i)=sum(sum(wJq.*(q1(1:Nq,:).*(Vq*rhs1) + q2(1:Nq,:).*(Vq*rhs2) + q3(1:Nq,:).*(Vq*rhs3))));
        end
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        rho = rho + rk4b(INTRK)*res1;
        m = m + rk4b(INTRK)*res2;
        E = E + rk4b(INTRK)*res3;        
    end                
    
    if mod(i,interval)==0 || i==Nsteps
        rhoq = Vq*rho;
        mq = Vq*m;
        Eq = Vq*E;
        q1 = Pq*V1(rhoq,mq,Eq);
        q2 = Pq*V2(rhoq,mq,Eq);
        q3 = Pq*V3(rhoq,mq,Eq);
        Usnap(:,sk) = [rho(:);m(:);E(:)];
        Vsnap(:,sk) = [q1(:);q2(:);q3(:)];
        sk = sk + 1;
    end
    
    if mod(i,5)==0 || i==Nsteps
        plot(xp,Vp*rho,'b-','linewidth',2)
        hold on
        plot(xp,Vp*m,'r-','linewidth',2)
        plot(xp,Vp*E,'k-','linewidth',2)
                
        % plot(xq,VqPq*V1(rhoq,mq,Eq),'b-','linewidth',2)
        % hold on
        % plot(xq,VqPq*V2(rhoq,mq,Eq),'r-','linewidth',2)
        % plot(xq,VqPq*V3(rhoq,mq,Eq),'k-','linewidth',2)
        
        title(sprintf('Time = %f, tstep %d out of %d',dt*i,i,Nsteps))
        ylim([-1 7])
        hold off
        drawnow
        
    end        
end

% semilogy(abs(rhstest))

%%

if opt==0
    % sine solution
    t = FinalTime;
    rhoex = @(x) (2 + sin(pi*(x - t)));
    uex = @(x) ones(size(x));
    pex = @(x) ones(size(x));
    mex = @(x) rhoex(x).*uex(x);
    Eex = @(x) pex(x)./(gamma-1) + .5*rhoex(x).*uex(x).^2;
    
    e1 = rhoex(xq) - Vq*rho;
    e2 = mex(xq) - Vq*m;
    e3 = Eex(xq) - Vq*E;
    L2err = sqrt(sum(sum(wJq.*(e1.^2+e2.^2+e3.^2))));
    L2err
end
%%

function [rhs1 rhs2 rhs3] = rhsEuler(rhoq,mq,Eq,rho,m,E,i,INTRK)

Globals1D
global rhoL rhoR mL mR EL ER gamma
global mapP mapB
global f1 f2 f3
global DNr Pq Lq
global Vf wq Vq

uq = mq./rhoq;
% pq = (gamma-1)*(Eq-.5*rhoq.*uq.^2);

rhoM = rhoq(end-1:end,:);
uM = uq(end-1:end,:);
EM = Eq(end-1:end,:);
mM = rhoM.*uM;
rhoP = reshape(rhoM(mapP),Nfp*Nfaces,K);
uP = reshape(uM(mapP),Nfp*Nfaces,K);
EP = reshape(EM(mapP),Nfp*Nfaces,K);
mP = rhoP.*uP;
pM = (gamma-1)*(EM-.5*rhoM.*uM.^2);

global useBC
if useBC==1   
    rhoP(mapB) = [rhoL rhoR];
    uP(mapB) = [mL/rhoL mR/rhoR];
    mP(mapB) = [mL mR];
    EP(mapB) = [EL ER];
elseif useBC==2    
    rhoP(mapB) = rhoM(mapB);
    EP(mapB) = EM(mapB);    
    uP(mapB) = -uM(mapB);
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
    
    FS = f1(rhox,rhoy,ux,uy,Ex,Ey);
    FS = rx(1,e)*sum(DNr.*FS,2);
    rhs1(:,e) = [Pq Lq]*FS;
    
    FS = f2(rhox,rhoy,ux,uy,Ex,Ey);
    FS = rx(1,e)*sum(DNr.*FS,2);
    rhs2(:,e) = [Pq Lq]*FS;
    
    FS = f3(rhox,rhoy,ux,uy,Ex,Ey);
    FS = rx(1,e)*sum(DNr.*FS,2);
    rhs3(:,e) = [Pq Lq]*FS; 
end

global tau

% local lax penalty
fluxopt = 1;
if fluxopt==1
    
    cvel = sqrt(gamma*pM./rhoM);
    lm   = (abs(uM) + cvel);
    Lfc  = max(lm(mapP),lm);
 
    d1 = Lfc.*((rhoP-rhoM));
    d2 = Lfc.*((mP-mM));
    d3 = Lfc.*((EP-EM));   
    
else
    
    global V1 V2 V3 avg pfun beta    
    du1 = (V1(rhoP,mP,EP)-V1(rhoM,mM,EM));
    du2 = (V2(rhoP,mP,EP)-V2(rhoM,mM,EM));
    du3 = (V3(rhoP,mP,EP)-V3(rhoM,mM,EM));
    u2avg = avg(uM.^2,uP.^2);
    uavg = avg(uM,uP);
    betaM = beta(rhoM,uM,EM);
    betaP = beta(rhoP,uP,EP);
    betalog = logmean(betaM,betaP);
    ubar = (2*uavg.^2-u2avg);
    pavg = avg(rhoM,rhoP)./(2*avg(betaM,betaP));
    rholog = logmean(rhoM,rhoP);
    h = gamma./(2*betalog*(gamma-1)) + .5*ubar;
    a = sqrt(gamma.*pavg./rholog);
    
    R11 = 1; R12 = 1; R13 = 1;
    R21 = uavg-a; R22 = uavg; R23 = uavg + a;
    R31 = h-uavg.*a; R32 = .5*ubar; R33 = h+uavg.*a;
    
    D11 = abs((uavg - a)).*rholog/(2*gamma);
    D22 = abs((uavg)).*rholog*(gamma-1)/gamma;
    D33 = abs((uavg + a)).*rholog/(2*gamma);
    
    r1 = D11.*(R11.*du1 + R21.*du2 + R31.*du3);
    r2 = D22.*(R12.*du1 + R22.*du2 + R32.*du3);
    r3 = D33.*(R13.*du1 + R23.*du2 + R33.*du3);
    
    d1 = R11.*r1 + R12.*r2 + R13.*r3;
    d2 = R21.*r1 + R22.*r2 + R23.*r3;
    d3 = R31.*r1 + R32.*r2 + R33.*r3;
end

f1f = -nx.*f1f + .5*tau*d1;
f2f = -nx.*f2f + .5*tau*d2;
f3f = -nx.*f3f + .5*tau*d3;
rhs1 = -2*rhs1 + Lq*(Fscale.*f1f);
rhs2 = -2*rhs2 + Lq*(Fscale.*f2f);
rhs3 = -2*rhs3 + Lq*(Fscale.*f3f);

end
