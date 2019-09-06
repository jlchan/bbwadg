clear

N = 9;
K = 16; % number of elements
CFL = .75; % should be O(1) - too large, soln blows up
FinalTime = 4;

VX = linspace(-1,1,K+1)';
h = VX(2)-VX(1);


%% nodal DG

[r w] = JacobiGQ(0,0,N);
V = Vandermonde1D(N,r);
Dr = GradVandermonde1D(N,r)/V;

[rq wq] = JacobiGQ(0,0,N);
Vq = Vandermonde1D(N,rq)/V;
M = Vq'*diag(wq)*Vq;

Vf = Vandermonde1D(N,[-1;1])/V;
L = M\(Vf');

Qr = M*Dr;
Qrskew = (Qr-Qr');

QNr = [Qrskew Vf'*diag([-1;1])
    -diag([-1;1])*Vf zeros(2)];
WN = diag([w;1;1]);
DNrskew = WN\QNr;
PN = M\[Vq' Vf']*WN;

J = h/2;

x = zeros(N+1,K);
for e = 1:K
    x(:,e) = h*(1+r)/2 + VX(e);
end

mapP = zeros(2,K);
for e = 2:K-1
    mapP(:,e) = [-2;1] + 2*e;
end
mapP(:,1) = [2*K; 3];
mapP(:,end) = [2*K-2; 1];

nx = repmat([-1;1],1,K);

%% make fv sub-cells

rf = [-1;cumsum(w)-1];
dx = h*repmat(w,1,K)/sum(w);

% map for L/R cell centers
mapL = reshape([(N+1)*K 1:(N+1)*K-1],N+1,K);
mapR = reshape([2:(N+1)*K 1],N+1,K);

mapPfv = zeros(2,(N+1)*K);
for e = 2:(N+1)*K-1
    mapPfv(:,e) = [-2;1] + 2*e;
end
mapPfv(:,1) = [2*(N+1)*K; 3];
mapPfv(:,end) = [2*(N+1)*K-2; 1];

nxfv = repmat([-1;1],1,(N+1)*K);

%%

tau = 1; % try tau b/w 0 and 1

% advection
fS = @(uL,uR) .5*(uL + uR);
fD = @(uL,uR) .5*(uR-uL);

% Burgers
fS = @(uL,uR) (uL.^2 + uL.*uR + uR.^2)/6;
fD = @(uL,uR) .5*abs(.5*(uL+uR)).*(uR-uL);

% interface flux
fSnfv = @(uL,uR,nx) fS(uL,uR).*nx - tau*fD(uL,uR); %.5*(uL.^b+uR.^b).*nx - tau*.5*(uR-uL).*abs(.5*(uR+uL)).^(b-1);

% initial guess
% u0 = @(x) -sin(pi*x);
u0 = @(x) exp(-50*x.^2); % gaussian initial condition
a = 100; u0 = @(x) 1e-15 + .5*(tanh(a*(x+.25))-tanh(a*(x-.25)));
% u0 = @(x) 1e-15 + (1-5*abs(x)).*(abs(x)<1/5);

% u0 = @(x) exp(-200*(x+.5).^2) + (1-10*abs(x)).*(abs(x)<1/10) + (abs(x-.5)<.1);
% u0 = @(x) 1 + sin(2*pi*x);

u = u0(x);
% plot(x,u,'o');return

ssprku = 1-[0 1/4 2/3];
ssprkb = [1 .25 2/3];

CN = (N+1)^2;
dt = CFL*min(dx(:))/(max(abs(u(:)))); % set time-step
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

figure(1)

unorm = zeros(Nsteps,1);
rhsDG = zeros(N+1,K);
for i = 1:Nsteps
    
    uu = u;
    for INTRK = 1:3
        
        % 2nd order FV
        dxL = .5*(dx + dx(mapL));
        dxR = .5*(dx(mapR) + dx);
        
        sL = (uu - uu(mapL))./(dxL); % backwards diff slopes
        sR = (uu(mapR) - uu)./(dxR); % forward diff slopes
        %         sC = (uu(mapR)-uu(mapL))./(dxL + dxR);
        %         slopes = minmod([2*sL(:)'; sC(:)'; 2*sR(:)'])'; % MUSCL
        slopes = minmod([sL(:)'; sR(:)'])'; % minmod
        uL = uu(:) - slopes .* dx(:)/2;
        uR = uu(:) + slopes .* dx(:)/2;
        uM = [uL';uR'];
        uP = uM(mapPfv);
        rhsfv = -reshape(sum(fSnfv(uM,uP,nxfv),1),N+1,K)./dx;
        
        % DG
        uM = Vf*uu;
        uP = uM(mapP);
        uN = [uu;uM];        
        flux = fSnfv(uM,uP,nx);
        for e = 1:K
            [ux uy] = meshgrid(uN(:,e));
            rhsDG(:,e) = -(PN*sum(DNrskew.*fS(ux,uy),2) + L*flux(:,e))./J;
        end
                
        uufv = ssprku(INTRK)*u + (1-ssprku(INTRK))*uu + ssprkb(INTRK)*dt*rhsfv;
        uuDG = ssprku(INTRK)*u + (1-ssprku(INTRK))*uu + ssprkb(INTRK)*dt*rhsDG;
                
        a = max(min(1,min(uufv)./(min(uufv)-min(uuDG))),0);
        a(min(u + dt*rhsDG) > 0) = 1;
        a = repmat(a,N+1,1);
%         a = 1;
        
        uu = uufv + a.*(uuDG-uufv);
        
    end
    u = uu;
    
    unorm(i) = sum(sum(dx.*u.^2));
    
    if i==Nsteps || (mod(i,50)==0)
        plot(x,u,'o--')
        title(sprintf('min u = %g\n',min(u(:))))
        axis([-1 1 -.5 1.5])
        drawnow
    end
end

% plot(x,u-u0(x),'o')

figure(2)
hold on
plot(x,u,'x--')
