clear
clear -global
Globals1D;

CFL = .75;
% CFL = .125/2.5;
% CFL = .125/4;
N = 4;
K1D = 8;
computeEigs = 0;
useSBP = 0;

FinalTime = 10;
global tau
tau = 0;

r = JacobiGL(0,0,N);

[rq wq] = JacobiGL(0,0,N); rq = rq*(1-1e-11);
[rq wq] = JacobiGL(0,0,N+3); rq = rq*(1-1e-11);
% [rq wq] = JacobiGQ(0,0,N+1);

Nq = length(rq);

[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

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

rp = linspace(-1,1,50); F = eye(N+1);
Vp = Vandermonde1D(N,rp)/V;
xp = Vandermonde1D(N,rp)/Vandermonde1D(N,JacobiGL(0,0,N))*x;

% filtering
global F
d = ones(N+1,1);
% d(1:3) = [1 .5 .25];
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
mapP(end) = 1; mapP(1) = Nfp*Nfaces*K;
mapP = reshape(mapP,Nfp*Nfaces,K);
% mapB = [1; Nfp*Nfaces*K];

global f1
avg = @(uL,uR) .5*(uL+uR);
f1 = @(uL,uR,aL,aR) avg(uL,uR).*avg(aL,aR); % D(au) + a*Du + u*Da
% f1 = @(uL,uR,aL,aR) uL.^2 + uL.*uR + uR.^2;
% f1 = @(uL,uR,aL,aR) aR.*uR; %.5*(uL.*aL + aR.*uL + aL.*uR + aR.*uR);

%%
x0 = 0;
uq = (2 + (abs(xq-x0) < .5));
uq = exp(-25*(xq+.5).^2);
uq = -sin(pi*xq);
u = Pq*(uq);

global aq
aq = [Vq;Vf]*Pq*(1 + exp(-10^2*(xq-.5).^2));
% aq = [Vq;Vf]*Pq*(1 + exp(-10^2*cos(pi*xq).^2));
% plot([xq;x([1 end],:)],aq,'o');return

%% get eigs

if computeEigs
    rxq = repmat(rx(1,:),size(Vq,2),1);
    nx = repmat([-1;1],1,K);    
    A = zeros(size(Vq,2)*K);
    uu = zeros(size(Vq,2),K);
    for i = 1:size(Vq,2)*K
        uu(i) = 1;
        uq = [Vq;Vf]*uu;
        [rhs] = rhsEuler(uq,i,1);
        A(:,i) = rhs(:);
        uu(i) = 0;
    end
    % S = kron(diag(J(1,:)),Vq'*diag(wq)*Vq)*A;
    [W D] = eig(A);
    lam = diag(D);
    plot(lam,'o')
    xlim([-.1 .1])
    title(sprintf('max real part = %g\n',max(real(lam))))
    return    
end

D = [Pq Lq]*DNr*[Vq; Vf];
daq = D*Pq*aq(1:Nq,:);
[~,id] = max(abs(daq(:)));


%%
res1 = 0;

CN = (N+1)^2/2;
CN = 1/min(wq);
h = 2/K1D;
dt = CFL * h/CN;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

global uL uR mL mR EL ER
global useBC; useBC = 1;


VqPq2 = [zeros(1,Nq+2);zeros(Nq,1) VqPq zeros(Nq,1);zeros(1,Nq+2)]; VqPq2(1) = 1; VqPq2(end) = 1;
xq2 = Vandermonde1D(N,[-1;rq;1])/V*x;
% plot(xq2,VqPq2*([Vf(1,:);Vq;Vf(2,:)]*u)); return

wJq = diag(wq)*((Vandermonde1D(N,rq)/V)*J);



%%

dt0 = dt;

figure(1)
for i = 1:Nsteps
    
    for INTRK = 1:5
        
        % interpolate to quadrature
        uq = [Vq;Vf]*u;
        
        if INTRK==1 && mod(i-1,5)==0
            uq1 = uq([Nq+1,1:Nq,Nq+2],:);
            clf
            hold on
            plot(xq2,uq1,'ro--','linewidth',2)
            
            wq0 = [0;wq;0];
            plot(.5*wq0'*xq2,.5*wq0'*uq1,'s','linewidth',2)
            plot(xq(id),xq(id)*0,'ro','markersize',32)
            axis([-1 1 -2 2])
            title(sprintf('timestep %d out of %d, RK step %d,time %f',i,Nsteps,INTRK,dt*i))
            
            if (i==1 && INTRK==1)
                pause
            else
                drawnow
            end
        end
        
        [rhs1] = rhsEuler(uq,i,INTRK);
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        u = u + rk4b(INTRK)*res1;
        
        %         energy(i) = sum(sum(wJq.*(aq(1:Nq,:).*(Vq*u).^2)));
        energy(i) = sum(sum(wJq.*((Vq*u).^2)));
        
    end
    
end

figure
semilogy(energy,'--')


function [rhs1] = rhsEuler(uq,i,INTRK)

Globals1D
global uL uR
global mapP mapB
global f1 f2 f3
global DNr Pq Lq
global Vf wq Vq
global aq

uM = uq(end-1:end,:);
uP = reshape(uM(mapP),Nfp*Nfaces,K);

aM = aq(end-1:end,:);
aP = reshape(aM(mapP),Nfp*Nfaces,K);

% compute fluxes
f1f = f1(uM,uP,aM,aP);

nhat = [-1; 1];
rhs1 = zeros(Np,K);
for e = 1:K
    [ux uy] = meshgrid(uq(:,e));
    [ax ay] = meshgrid(aq(:,e));
    
    FS = f1(ux,uy,ax,ay);
    fu = f1(uM(:,e),uM(:,e),aM(:,e),aM(:,e));
    FS = rx(1,e)*sum(DNr.*FS,2);
    fS = nhat.*(f1f(:,e) - fu);
    rhs1(:,e) = [Pq Lq]*FS + Lq*(.5*Fscale(1,e)*fS);
end

global tau
% local lax penalty
Lfc = max(abs(aM),abs(aP));
d1 = Lfc.*((uP-uM));

rhs1 = -2*rhs1 + .5*tau*Lq*(Fscale.*d1);


end