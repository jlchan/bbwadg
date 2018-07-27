clear -global
Globals1D;

projectV = 1;
CFL = .125;
% CFL = .125/2.5;
N = 4;
K1D = 16;

useSBP = 0;

FinalTime = 1.8;
FinalTime = .2;
opt = 2;

global tau
tau = 1;

% r = JacobiGL(0,0,N);
r = JacobiGL(0,0,N);

[rq wq] = JacobiGL(0,0,N); rq = rq*(1-1e-11);
% [rq wq] = JacobiGL(0,0,N+1); rq = rq*(1-1e-11);
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

global VqPq Ef Dq Vq Pq Dfq Lq
global Vf wq Vq
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
DNr = [Drq-.5*Vq*Lq*diag(nrJ)*VfPq .5*VqLq*diag(nrJ); -.5*diag(nrJ)*VfPq .5*diag(nrJ)*eye(2)];
W = diag([wq;1;1]);

Ef = zeros(2,length(rq)); Ef(1) = 1; Ef(end) = 1;

global Vp xp
rp = linspace(-1,1,25); F = eye(N+1);
Vp = Vandermonde1D(N,rp)/V;
xp = Vandermonde1D(N,rp)/Vandermonde1D(N,JacobiGL(0,0,N))*x;

% filtering
global F
d = ones(N+1,1);
d(end) = 0;
F = V*(diag(d)/V);

%% switch to SBP

if useSBP
    % for SBP operator (with boundary nodes)
    Np = length(rq);
    Pq = eye(length(rq));
    Vq = eye(length(rq));
    Vf = zeros(2,length(rq)); Vf(1,1) = 1; Vf(end,end) = 1;
    Lq = diag(1./wq)*Vf';
    VqPq = Vq*Pq;
    VfPq = Vf*Pq;
    VqLq = Vq*Lq;
end

% check SBP property
Iq = eye(length(wq));
Ef = zeros(2,length(rq)); Ef(1) = 1; Ef(end) = 1;
Qsbp = [Iq; Ef]'*W*DNr*[Iq;Ef];
% Qsbp+Qsbp'

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
    Eq = pq/(gamma-1) + .5*rhoq.*uq.^2;
    mq = rhoq.*uq;
    
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
    
    rhoRex = @(x) 1 + .2*sin(5*x);
    rhoL = 3.857143; rhoR = rhoRex(x(end));
    uL   = 2.629369; uR = 0;
    pL   = 10.3333;  pR = 1;
    
    %     rhoRex = @(x) 2 + .2*sin(5*x);
    %     rhoL = 4; rhoR = rhoRex(x(end));
    %     uL   = 2; uR = 0;
    %     pL   = 4; pR = 2;
    
    %     useBC = 0;
    %     mapP(1) = mapM(end); mapP(end) = mapM(1);
    %     vmapP(1) = vmapM(end); vmapP(end) = vmapM(1);
    
    rhoq = rhoL*(xq < -4) + (rhoRex(xq)).*(xq >= -4);
    uq = uL*(xq < -4);
    pq = pL*(xq < -4) + pR*(xq >= -4);
    
    mL = rhoL*uL; mR = rhoR*uR;
    EL = pL/(gamma-1) + .5*rhoL*uL.^2;
    ER = pR/(gamma-1) + .5*rhoR*uR.^2;
    
elseif opt==4 % modified sod
    
    rhoL = 1; rhoR = .125;
    uL = .75; uR = 0;
    pL = 1; pR = .1;
    
    rhoq = rhoL*(xq < .3) + rhoR*(xq > .3);
    uq = uL*(xq < .3);
    pq = pL*(xq < .3) + pR*(xq > .3);
    
    mL = rhoL*uL; mR = rhoR*uR;
    EL = pL/(gamma-1) + .5*rhoL*uL.^2;
    ER = pR/(gamma-1) + .5*rhoR*uR.^2;
    
end

rho = Pq*(rhoq);
u = Pq*(uq);
p = Pq*(pq);
E = Pq*(pq/(gamma-1) + .5*rhoq.*uq.^2);
m = Pq*(rhoq.*uq);

VqPq2 = [zeros(1,Nq+2);zeros(Nq,1) VqPq zeros(Nq,1);zeros(1,Nq+2)]; VqPq2(1) = 1; VqPq2(end) = 1;
xq2 = Vandermonde1D(N,[-1;rq;1])/V*x;
% plot(xq2,VqPq2*([Vf(1,:);Vq;Vf(2,:)]*rho)); return

wJq = diag(wq)*((Vandermonde1D(N,rq)/V)*J);

%% limiter stuff

global Sfun SU lim Sref
fids = size(Vq,1)-1:size(Vq,1);
VV = @(U) [V1(U(:,1),U(:,2),U(:,3)),V2(U(:,1),U(:,2),U(:,3)),V3(U(:,1),U(:,2),U(:,3))];
UU = @(V) [U1(V(:,1),V(:,2),V(:,3)),U2(V(:,1),V(:,2),V(:,3)),U3(V(:,1),V(:,2),V(:,3))];

rhoq = Vq*rho; mq = Vq*m; Eq = Vq*E;
Sref = max(max(-rhoq.*s(rhoq,mq,Eq)));
% s = @(rho,m,E) log((gamma-1).*rhoe(rho,m,E)./(rho.^gamma));
Sfun = @(rho,m,E) -rho.*s(rho,m,E) - Sref - 1;
lim = @(u,theta) .5*wq'*u + theta.*(u - .5*wq'*u);
SU = @(U) Sfun(U(:,1),U(:,2),U(:,3));

% plot(xp,Sfun(Vp*rho,Vp*m,Vp*E),'o')
% return

%% adaptive timestepping

if 1
    U0 = [rho(:);m(:);E(:)];
    
    options = odeset('RelTol',1e-3,'AbsTol',1e-6);
    [t,Uf] = ode45(@(t,U) rhs45(t,U),[0,FinalTime],U0,options);
    %     [t,Uf] = ode45(@(t,U) rhs45(t,U),[0,FinalTime],U0);
    
    U = Uf(end,:)';
    ids = 1:Np*K;
    rho = reshape(U(ids),Np,K);
    m = reshape(U(ids + Np*K),Np,K);
    E = reshape(U(ids + 2*Np*K),Np,K);
    
    rhoq45 = Vq*rho;
    mq45 = Vq*m;
    Eq45 = Vq*E;
    
    rhop45 = Vp*rho;
    mp45 = Vp*m;
    Ep45 = Vp*E;
    up45 = mp45./rhop45;
    pp45 = (gamma-1)*(Ep45-.5*rhop45.*up45.^2);
    plot(xp,rhop45,'b-','linewidth',2)
    hold on
    plot(xp,up45,'r-.','linewidth',2)
    plot(xp,pp45,'k--','linewidth',2)
    
    keyboard
    return
end
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
h = 2/K1D;
dt = CFL * h/CN;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

dt0 = dt;
time = 0;
tstep = 1;

figure(1)
while time < FinalTime
    i = tstep;
    for INTRK = 1:5
        
        % interpolate to quadrature
        rhoq = Vq*rho;
        mq = Vq*m;
        Eq = Vq*E;
        
        % project entropy variables
        q1 = [Vq;Vf]*Pq*V1(rhoq,mq,Eq);
        q2 = [Vq;Vf]*Pq*V2(rhoq,mq,Eq);
        q3 = [Vq;Vf]*Pq*V3(rhoq,mq,Eq);
        
        rhoq = U1(q1,q2,q3);
        mq   = U2(q1,q2,q3);
        Eq   = U3(q1,q2,q3);
        
        aa = 5;
        ratio = max(Eq,[],1)./max([Vq;Vf]*E,[],1);
        %         ratio = max(rhoq,[],1)./max([Vq;Vf]*rho,[],1);
        
        %         aa = 1.75;
        %         S0 = Sfun([Vq;Vf]*rho,[Vq;Vf]*m,[Vq;Vf]*E);
        % %         ratio = max(abs(Sfun(rhoq,mq,Eq)-S0)./Sref);
        % ratio = min(Sfun(rhoq,mq,Eq),[],1)./min(S0,[],1);
        
        %         if max(ratio) > aa
        %             xq2 = [Vq;Vf]*x;
        %             figure
        %             plot(xq2,Sfun(rhoq,mq,Eq),'o')
        %             hold on
        %
        %             % project entropy variables
        %             rhoq = Vq*rho;
        %             mq = Vq*m;
        %             Eq = Vq*E;
        %             theta = .5;
        %             q1 = [Vq;Vf]*Pq*V1(lim(rhoq,theta),lim(mq,theta),lim(Eq,theta));
        %             q2 = [Vq;Vf]*Pq*V2(lim(rhoq,theta),lim(mq,theta),lim(Eq,theta));
        %             q3 = [Vq;Vf]*Pq*V3(lim(rhoq,theta),lim(mq,theta),lim(Eq,theta));
        %
        %             rhoq = U1(q1,q2,q3);
        %             mq   = U2(q1,q2,q3);
        %             Eq   = U3(q1,q2,q3);
        %             plot(xq2,Sfun(rhoq,mq,Eq),'x')
        %             keyboard
        %
        %             e = find(ratio > aa,1,'first');
        %             U = [rho(:,e) m(:,e) E(:,e)];
        %             [theta fval] = bisect_lim(U,aa);
        %         end
        
        if 0 % limiter
            %             fprintf('limiting on tstep %d\n',i);
            elist = find(ratio>aa);
            elistK{i} = elist;
            for ee = 1:length(elist)
                
                e = elist(ee);
                U = [rho(:,e) m(:,e) E(:,e)];
                
                [theta fval] = bisect_lim(U,aa);
                U = lim(Vq*U,theta);
                
                rho(:,e) = Pq*U(:,1);
                m(:,e) = Pq*U(:,2);
                E(:,e) = Pq*U(:,3);
                
                % re-project entropy variables
                q1(:,e) = [Vq;Vf]*Pq*V1(U(:,1),U(:,2),U(:,3));
                q2(:,e) = [Vq;Vf]*Pq*V2(U(:,1),U(:,2),U(:,3));
                q3(:,e) = [Vq;Vf]*Pq*V3(U(:,1),U(:,2),U(:,3));
            end
            rhoq = U1(q1,q2,q3);
            mq   = U2(q1,q2,q3);
            Eq   = U3(q1,q2,q3);
        end
        
        rhoq1 = rhoq([Nq+1,1:Nq,Nq+2],:);
        rhoq2 = [Vf(1,:);Vq;Vf(2,:)]*rho;
        mq1 = mq([Nq+1,1:Nq,Nq+2],:);
        mq2 = [Vf(1,:);Vq;Vf(2,:)]*m;
        Eq1 = Eq([Nq+1,1:Nq,Nq+2],:);
        Eq2 = [Vf(1,:);Vq;Vf(2,:)]*E;
        
        if INTRK==5 && mod(i-1,5)==0
            clf
            %             subplot(2,1,1)
            hold on
            plot(xq2,rhoq1,'ro--','linewidth',2)
            plot(xq2,rhoq2,'b-','linewidth',2)
            
            plot(xq2,Eq1,'g^--','linewidth',2)
            plot(xq2,Eq2,'k-','linewidth',2)
            
            wq0 = [0;wq;0];
            plot(.5*wq0'*xq2,.5*wq0'*rhoq1,'s','linewidth',2)
            plot(.5*wq0'*xq2,.5*wq0'*rhoq2,'x','linewidth',2)
            plot(.5*wq0'*xq2,.5*wq0'*Eq1,'s','linewidth',2)
            plot(.5*wq0'*xq2,.5*wq0'*Eq2,'x','linewidth',2)
            
            %             subplot(2,1,2)
            %             hold on
            %             plot(xq2,V1(rhoq2,mq2,Eq2),'ro--','linewidth',2)
            %             plot(xq2,[Vf(1,:);Vq;Vf(2,:)]*Pq*V1(rhoq(1:Nq,:),mq(1:Nq,:),Eq(1:Nq,:)),'b-','linewidth',2)
            %             plot(xq2,V3(rhoq2,mq2,Eq2),'g^--','linewidth',2)
            %             plot(xq2,[Vf(1,:);Vq;Vf(2,:)]*Pq*V3(rhoq(1:Nq,:),mq(1:Nq,:),Eq(1:Nq,:)),'k-','linewidth',2)
            
            title(sprintf('timestep %d out of %d, RK step %d,time %f',i,Nsteps,INTRK,dt*i))
            %                         axis([-.5 .5 1e-11 .25])
            % pause
            
            if (i==1 && INTRK==1)
                pause
            else
                drawnow
            end
            
            
        end
        
        [rhs1 rhs2 rhs3] = rhsEuler(rhoq,mq,Eq,rho,m,E,i,INTRK);
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        rho = rho + rk4b(INTRK)*res1;
        m = m + rk4b(INTRK)*res2;
        E = E + rk4b(INTRK)*res3;
        
    end
    
    rhoq = Vq*rho;
    mq = Vq*m;
    Eq = Vq*E;
    Sq = Sfun(rhoq,mq,Eq);
    Si = sum(sum(wJq.*Sq));
    Svec(i) = Si;
    
    if 0 % mod(i,5)==0 || i==Nsteps
        uq = mq./rhoq;
        pq = (gamma-1)*(Eq-.5*rhoq.*uq.^2);
        plot(xp,Vp*Pq*rhoq,'b-','linewidth',2)
        hold on
        %         plot(xq,VqPq*uq,'r-','linewidth',2)
        %         plot(xq,VqPq*pq,'k-','linewidth',2)
        
        plot(.5*wq'*xq,.5*wq'*rhoq,'s','linewidth',2)
        plot(.5*wq'*xq,.5*wq'*uq,'s','linewidth',2)
        plot(.5*wq'*xq,.5*wq'*pq,'s','linewidth',2)
        
        % plot(xq,VqPq*V1(rhoq,mq,Eq),'b-','linewidth',2)
        % hold on
        % plot(xq,VqPq*V2(rhoq,mq,Eq),'r-','linewidth',2)
        % plot(xq,VqPq*V3(rhoq,mq,Eq),'k-','linewidth',2)
        
        title(sprintf('Time = %f, entropy = %f',dt*i,Si))
        %         axis([-5 5 -1 7])
        hold off
        drawnow
        
    end
    
    time = time + dt;
    tstep = tstep + 1;
end

figure(2)
semilogy(dt*(1:tstep-1),abs(Svec(1:tstep-1)-Svec(1)),'--')
hold on


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

function rhs = rhs45(t,U)

Globals1D
global VqPq tf Dq Vq Vf Pq Dfq
global Vp xp
global V1 V2 V3 U1 U2 U3
ids = 1:Np*K;
rho = reshape(U(ids),Np,K);
m = reshape(U(ids + Np*K),Np,K);
E = reshape(U(ids + 2*Np*K),Np,K);

% interpolate to quadrature
rhoq = Vq*rho;
mq = Vq*m;
Eq = Vq*E;

% project entropy variables
q1 = [Vq;Vf]*Pq*V1(rhoq,mq,Eq);
q2 = [Vq;Vf]*Pq*V2(rhoq,mq,Eq);
q3 = [Vq;Vf]*Pq*V3(rhoq,mq,Eq);
rhoq = U1(q1,q2,q3);
mq   = U2(q1,q2,q3);
Eq   = U3(q1,q2,q3);

if any(isnan([rhoq(:);mq(:);Eq(:)])) || any(norm(imag([rhoq(:);mq(:);Eq(:)]),'fro')>1e-10)
    keyboard
end

[rhs1 rhs2 rhs3] = rhsEuler(rhoq,mq,Eq);
rhs = [rhs1(:);rhs2(:);rhs3(:)];

end

function [rhs1 rhs2 rhs3] = rhsEuler(rhoq,mq,Eq,rho,m,E,i,INTRK)

Globals1D
global rhoL rhoR mL mR EL ER gamma
global mapP mapB
global f1 f2 f3
global DNr Pq Lq
global Vf wq Vq

uq = mq./rhoq;
pq = (gamma-1)*(Eq-.5*rhoq.*uq.^2);

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

nhat = [-1; 1];
rhs1 = zeros(Np,K);
rhs2 = zeros(Np,K);
rhs3 = zeros(Np,K);
for e = 1:K
    [rhox rhoy] = meshgrid(rhoq(:,e));
    [ux uy] = meshgrid(uq(:,e));
    [Ex Ey] = meshgrid(Eq(:,e));
    
    FS = f1(rhox,rhoy,ux,uy,Ex,Ey);
    fu = f1(rhoM(:,e),rhoM(:,e),uM(:,e),uM(:,e),EM(:,e),EM(:,e));
    FS = rx(1,e)*sum(DNr.*FS,2);
    fS = nhat.*(f1f(:,e) - fu);
    rhs1(:,e) = [Pq Lq]*FS + Lq*(.5*Fscale(1,e)*fS);
    
    FS = f2(rhox,rhoy,ux,uy,Ex,Ey);
    fu = f2(rhoM(:,e),rhoM(:,e),uM(:,e),uM(:,e),EM(:,e),EM(:,e));
    FS = rx(1,e)*sum(DNr.*FS,2);
    fS = nhat.*(f2f(:,e) - fu);
    rhs2(:,e) = [Pq Lq]*FS + Lq*(.5*Fscale(1,e)*fS);
    
    FS = f3(rhox,rhoy,ux,uy,Ex,Ey);
    fu = f3(rhoM(:,e),rhoM(:,e),uM(:,e),uM(:,e),EM(:,e),EM(:,e));
    FS = rx(1,e)*sum(DNr.*FS,2);
    fS = nhat.*(f3f(:,e) - fu);
    rhs3(:,e) = [Pq Lq]*FS + Lq*(.5*Fscale(1,e)*fS);
    
end

global tau

% local lax penalty
fluxopt = 1;
if fluxopt==1
    cvel = sqrt(gamma*pM./rhoM);
    lm   = (abs(uM) + cvel);
    Lfc  = max(lm(mapP),lm);
%          lm4 = abs(uM).^4 + cvel.^4;
%         Lfc = (.5*(lm4+lm4(mapP))).^(1/4);
    
    d1 = Lfc.*((rhoP-rhoM));
    d2 = Lfc.*((mP-mM));
    d3 = Lfc.*((EP-EM));
    
elseif fluxopt==2
    
    global V1 V2 V3
    cvel = sqrt(gamma*pM./rhoM);
    lm   = (abs(uM) + cvel);
        Lfc  = max(lm(mapP),lm);
    lm4 = abs(uM).^4 + cvel.^4;
    Lfc = (.5*(lm4+lm4(mapP))).^(1/4);
    
    v1P = V1(rhoP,mP,EP);    v1M = V1(rhoM,mM,EM);
    v2P = V2(rhoP,mP,EP);    v2M = V2(rhoM,mM,EM);
    v3P = V3(rhoP,mP,EP);    v3M = V3(rhoM,mM,EM);
    
    d1 = Lfc.*((v1P-v1M));
    d2 = Lfc.*((v2P-v2M));
    d3 = Lfc.*((v3P-v3M));
    
else
    
    global V1 V2 V3
    du1 = (V1(rhoP,mP,EP)-V1(rhoM,mM,EM));
    du2 = (V2(rhoP,mP,EP)-V2(rhoM,mM,EM));
    du3 = (V3(rhoP,mP,EP)-V3(rhoM,mM,EM));
    
    
    global avg pfun beta
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

rhs1 = -2*rhs1 + .5*tau*Lq*(Fscale.*d1);
rhs2 = -2*rhs2 + .5*tau*Lq*(Fscale.*d2);
rhs3 = -2*rhs3 + .5*tau*Lq*(Fscale.*d3);

end


function [theta ff] = bisect_lim(U,aa)

Globals1D
global U1 U2 U3 V1 V2 V3 SU Sfun
global Vf Pq Vq wq lim rp Vp xp

a = 0; b = 1;
ff0 = bisect_lim_func(U,a,aa);
ff1 = bisect_lim_func(U,b,aa);

if ff0*ff1 > 0
    % no root or blowup!
    figure
    for theta = 0:.01:1
        plot(theta,bisect_lim_func(U,theta,aa),'o')
        hold on
    end
    
    VU = @(U) [V1(U(:,1),U(:,2),U(:,3)),V2(U(:,1),U(:,2),U(:,3)),V3(U(:,1),U(:,2),U(:,3))];
    UV = @(V) [U1(V(:,1),V(:,2),V(:,3)),U2(V(:,1),V(:,2),V(:,3)),U3(V(:,1),V(:,2),V(:,3))];
    
    keyboard
end

for i = 1:10
    c = .5*(a+b);
    ff = bisect_lim_func(U,c,aa);
    if ff0*ff < 0
        b = c;
        ff1 = ff;
    elseif ff1*ff < 0
        a = c;
        ff0 = ff;
    end
end
theta = c;

end

function f = bisect_lim_func(U,theta,aa)

global U1 U2 U3 V1 V2 V3 SU Sfun
global Vf Pq Vq wq lim
global Sref

Uqf = [Vq;Vf]*U;
rholim = lim(Vq*U(:,1),theta);
mlim = lim(Vq*U(:,2),theta);
Elim = lim(Vq*U(:,3),theta);

q1 = [Vq;Vf]*Pq*V1(rholim,mlim,Elim);
q2 = [Vq;Vf]*Pq*V2(rholim,mlim,Elim);
q3 = [Vq;Vf]*Pq*V3(rholim,mlim,Elim);

f = max(U3(q1,q2,q3)) - aa*max(Vf*U(:,3));

% % spikes = negative
% % negcut = @(x) x.*(x<1e-10);
% S2 = Sfun(U1(q1,q2,q3),U2(q1,q2,q3),U3(q1,q2,q3));
% S1 = Sfun(Uqf(:,1),Uqf(:,2),Uqf(:,3));
% f = min(S2)-aa*min(S1);


end