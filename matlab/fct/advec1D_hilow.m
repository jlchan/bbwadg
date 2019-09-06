clear
clear -global
Globals1D;

CFL = .5;

N = 15;
K1D = 2;

FinalTime = 2;

r = JacobiGL(0,0,N);
[rq wq] = JacobiGL(0,0,N); rq = rq*(1-1e-11);
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

global rxJ sJ wq
rxJ = rx.*J;
sJ = ones(size(Fscale));

V = Vandermonde1D(N,r);
Vq = Vandermonde1D(N,rq)/V;
Dr = GradVandermonde1D(N,r)/V;
Vf = Vandermonde1D(N,[-1 1])/V;

M = Vq'*diag(wq)*Vq;
LIFT = M\Vf';

rp = linspace(-1,1,50);
Vp = Vandermonde1D(N,rp)/V;
xp = Vandermonde1D(N,rp)/V*x;

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


%% extract matricees

S = zeros((N+1)*K);
B = zeros((N+1)*K);
u = zeros(N+1,K);
for i = 1:(N+1)*K
    u(i) = 1;
    [rhs rS rB] = advec1D(u,1);
    u(i) = 0;
    S(:,i) = rS(:);
    B(:,i) = rB(:);
end
S(abs(S)<1e-8) =0;
B(abs(B)<1e-8) =0;

S = sparse(S);
B = sparse(B);

e = ones((N+1)*K,1);
SL = diag(e(2:end),1) - diag(e(2:end),-1);
SL(1,end) = -1;
SL(end,1) = 1;
SL = SL/2;
BL = abs(SL);

%%

u = .5+sin(pi*x);
u = abs(x)<.25;
% u = exp(-100*x.^2);

% rc = [-.5;.5];
% xc = Vandermonde1D(N,rc)/V * x;
% u = sin(pi*xc);

CN = 1/min(wq);
h = 2/K1D;
dt = CFL * h/CN;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;
wJq = diag(wq)*(Vq*J);

ssprk = [0 3/4 1/3];

u = u(:);
wJq = wJq(:);
x = x(:);

figure(1)
uH = u;
uL = u;
for i = 1:Nsteps
    
    uold = uH;
    for INTRK = 1:3
        
        [ux uy] = meshgrid(uH);
        FH = -(S.*(ux+uy) - B.*(ux-uy));
        uH = uH + dt*sum(FH,2)./wJq;
        
        [ux uy] = meshgrid(uL);
        FL = -(SL.*(ux+uy) - BL.*(ux-uy));
        uL = uL + dt*sum(FL,2)./wJq;
        
        uH = ssprk(INTRK)*uold + (1-ssprk(INTRK))*uH; 
        uL = ssprk(INTRK)*uold + (1-ssprk(INTRK))*uL; 
    end
    
    unorm(i) = sum(wJq.*uH.^2);
    umean(:,i) = sum(reshape(wJq.*uH,N+1,K),1);
    uLmean(:,i) = sum(reshape(wJq.*uL,N+1,K),1);
    
    if mod(i,10)==0
        %plot(xc(:),u,'o-')
        
        plot(x,uH,'.-')
        hold on
        plot(x,uL,'o-')
        hold off
        title(sprintf('min u = %g, max u = %g, mean(u) = %g, unorm(u) = %g\n',min(u),max(u),umean(i),unorm(i)))
        axis([-1,1,-1, 2])
        drawnow
    end
end

%%

figure(2)
% semilogy(unorm,'o--')
% semilogy(sum(abs(umean),1),'x--');hold on; semilogy(sum(abs(uLmean),1),'o--')
semilogy(abs(umean'),'x--'); hold on; semilogy(abs(uLmean'),'o--')
ylim([1e-7 1])
% hold on
% plot(x,u,'.-')
% title(sprintf('min u = %g, max u = %g, mean(u) = %g, unorm(u) = %g\n',min(u),max(u),sum(wJq.*u),unorm(i)))
% axis([-1,1,-1, 2])
%%


function [rhs rS rB] = advec1D(u,dt)

Globals1D

global mapP mapB
global rxJ sJ wq

uM = u([1 N+1],:);
uP = reshape(uM(mapP),Nfp*Nfaces,K);
flux = .5*(uP-uM).*nx;

% compute fluxes
rhs = rxJ.*(Dr*u) + LIFT*(sJ.*flux);
rS = diag(wq)*rhs;
rB = diag(wq)*LIFT*(sJ.*(uP-uM));

rhs = rhs - LIFT*(sJ.*(uP-uM));
rhs = -rhs./J;

end