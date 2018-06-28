clear -global
Globals1D;

projectV = 1;
CFL = .5;
% CFL = .125/2.5;
% CFL = .125/4;
N = 4;
K1D = 16;
FinalTime = 4;

useSBP = 0;

global tau
tau = 1;

r = JacobiGL(0,0,N);
% r = JacobiGQ(0,0,N);

% [rq wq] = JacobiGL(0,0,N); rq = rq*(1-1e-11);
% [rq wq] = JacobiGL(0,0,N+4); rq = rq*(1-1e-9);
[rq wq] = JacobiGL(0,0,N);

% % % include boundary nodes for extraction
% rq = [-1*.999999999999;rq;1*.999999999999];
% wq = [0;wq;0];
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
mapP = reshape(mapP,Nfp*Nfaces,K);
mapB = [1; Nfp*Nfaces*K];

%% fluxes

g = 1;
V1 = @(h,u) g*h - .5*u.^2;
V2 = @(h,u) u;

hfun = @(v1,v2) (v1+.5*v2.^2)/g;
ufun = @(v1,v2) v2;

global f1 f2
f1 = @(hL,hR,uL,uR) (hL.*uL + hR.*uR)/2;
f2 = @(hL,hR,uL,uR) .5*(hL.*uL + hR.*uR).*.5.*(uL+uR) + .5*g.*(hL.*hR);

Sfun = @(h,u) .5*(h.*u.^2 + g*h.^2);
%% solution

hex  = @(x) 2+cos(pi*x);
hex = @(x) 4 + exp(-25*x.^2);
hvex = @(x) 0*sin(pi*x);

% a = .25;
% hex  = @(x) a + 1*(x > 0);
% hvex = @(x) zeros(size(x));

h = Pq*hex(xq);
hu = Pq*hvex(xq);

%%
res1 = 0;
res2 = 0;

CN = (N+1)^2/2;
hh = 2/K1D;
dt = CFL * hh/CN;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

global rhoL rhoR mL mR EL ER
global useBC; 
% periodic
useBC = 0;
mapP(1) = mapM(end); mapP(end) = mapM(1);
vmapP(1) = vmapM(end); vmapP(end) = vmapM(1);

wJq = diag(wq)*((Vandermonde1D(N,rq)/V)*J);
VqPq2 = [zeros(1,Nq+2);zeros(Nq,1) VqPq zeros(Nq,1);zeros(1,Nq+2)]; VqPq2(1) = 1; VqPq2(end) = 1;
xq2 = Vandermonde1D(N,[-1;rq;1])/V*x;
P = zeros(Nq+2); p = [Nq+1, 1:Nq, Nq+2];
P(:,p) = eye(Nq+2);

figure(1)
for i = 1:Nsteps
    
    for INTRK = 1:5               
        
        % interpolate to quadrature
        hq = Vq*h; 
        huq = Vq*hu;                
        uq = huq./hq;
        
        % project entropy variables
        v1 = [Vq;Vf]*Pq*V1(hq,uq);
        v2 = [Vq;Vf]*Pq*V2(hq,uq);                         
        
        if projectV
            hq = hfun(v1,v2);
            uq = ufun(v1,v2);            
        end    
        
%         if mod(i,25)==0
%             clf
%             subplot(2,1,1)
%             hold on
%             plot(xp,Vp*h,'-','linewidth',2);
%             plot(xq2,P*hfun(v1,v2),'x--','linewidth',2);
%             plot(xp,(Vp*hu)./(Vp*h),'-','linewidth',2);
%             plot(xq2,P*ufun(v1,v2),'x--','linewidth',2);
%             
%             subplot(2,1,2)
%             semilogy(xq,1e-10+abs(Sfun(Vq*h,(Vq*hu)./(Vq*h))-Sfun(hq(1:end-2,:),uq(1:end-2,:))),'-','linewidth',2)
%             hold on
%             drawnow
%         end
        
        [rhs1 rhs2] = rhsSWE(hq,uq,i,INTRK);        
           
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        h = h + rk4b(INTRK)*res1;
        hu = hu + rk4b(INTRK)*res2;
        
    end
        
    hq = Vq*h;
    uq = (Vq*hu)./hq;
    S(i) = sum(sum(wJq.*(Sfun(hq,uq))));
                    
    if mod(i,5)==0 || i==Nsteps

        plot(xq,hq,'b-','linewidth',2)
        hold on
        plot(xq,uq,'r-','linewidth',2)
        
        title(sprintf('Time = %f',dt*i))
        %         axis([-5 5 -1 7])
        hold off
        drawnow
    end
end
semilogy(dt*(1:Nsteps),S,'--','linewidth',2)
axis([0 FinalTime mean(S)-1 mean(S)+1])

function [rhs1 rhs2] = rhsSWE(hq,uq,i,INTRK)

Globals1D
global rhoL rhoR mL mR EL ER gamma
global mapP mapB
global f1 f2 
global DNr Pq Lq 
global Vf wq Vq


hM = hq(end-1:end,:);
uM = uq(end-1:end,:); 
hP = reshape(hM(mapP),Nfp*Nfaces,K);
uP = reshape(uM(mapP),Nfp*Nfaces,K);

global useBC
if useBC
    hM(mapB) = [hL hR];
    uM(mapB) = [uL uR];    
end

% compute fluxes
f1f = f1(hM,hP,uM,uP);
f2f = f2(hM,hP,uM,uP);

nhat = [-1; 1];
rhs1 = zeros(Np,K);
rhs2 = zeros(Np,K);
for e = 1:K
    [hx hy] = meshgrid(hq(:,e));
    [ux uy] = meshgrid(uq(:,e));
        
    FS = f1(hx,hy,ux,uy);    
    fu = f1(hM(:,e),hM(:,e),uM(:,e),uM(:,e));    
    FS = rx(1,e)*sum(DNr.*FS,2);    
    fS = nhat.*(f1f(:,e) - fu);    
    rhs1(:,e) = [Pq Lq]*FS + Lq*(.5*Fscale(1,e)*fS);
        
    FS = f2(hx,hy,ux,uy);    
    fu = f2(hM(:,e),hM(:,e),uM(:,e),uM(:,e));    
    FS = rx(1,e)*sum(DNr.*FS,2);    
    fS = nhat.*(f2f(:,e) - fu); 
    rhs2(:,e) = [Pq Lq]*FS + Lq*(.5*Fscale(1,e)*fS);            
end

global tau

% local lax penalty
lM = abs(uM) + sqrt(hM);
Lfc = max(lM,lM(mapP));
d1 = Lfc.*((hP-hM));
d2 = Lfc.*((hP.*uP-hM.*uM));

rhs1 = -2*rhs1 + .5*tau*Lq*(Fscale.*d1);
rhs2 = -2*rhs2 + .5*tau*Lq*(Fscale.*d2);


end