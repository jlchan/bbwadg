% SBP spectral method

clear
N = 25;
x = linspace(-1,1,2*N+2)'; % interpolatory
x = x(1:end-1);
dx = x(2)-x(1);

sk = 1;
D = zeros(2*N+1);
for k = -N:N
    lamk = 1i*k*pi;
    Vq(:,sk) = exp(lamk*x)/sqrt(2);        
    D(sk,sk) = lamk;
    sk = sk + 1;
end
% plot(x,Vq)

M = real(dx*(Vq'*Vq));
W = dx*eye(2*N+1);
Pq = M\(Vq'*dx);
D = real(Vq*D*Pq); % imag part = 0

K = (D'*W*D); %/dx;
Q = W*D;

%% make 2D

Np = (2*N+1)^2;

[x y] = meshgrid(x);

%%

% Burgers
fS = @(uL,uR) (uL.^2 + uL.*uR + uR.^2)/6;
fD = @(uL,uR) max(abs(uL),abs(uR)).*(uR-uL);

%%

FinalTime = 2;

rk4a = [            0.0 ...
    -567301805773.0/1357537059087.0 ...
    -2404267990393.0/2016746695238.0 ...
    -3550918686646.0/2091501179385.0  ...
    -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
    5161836677717.0/13612068292357.0 ...
    1720146321549.0/2090206949498.0  ...
    3134564353537.0/4481467310338.0  ...
    2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
    1432997174477.0/9575080441755.0 ...
    2526269341429.0/6820363962896.0 ...
    2006345519317.0/3224310063776.0 ...
    2802321613138.0/2924317926251.0];

u = -sin(pi*x).*sin(pi*y);

dt = .5*dx;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;
res = zeros(size(u));
rhs = zeros(size(u));

interval = 1;
Usnap = zeros(Np,ceil(Nsteps/interval)+1);
Usnap(:,1) = u(:);
sk = 2;

for i = 1:Nsteps
    for INTRK = 1:5
        
        % x coordinate
        for ii = 1:2*N+1
            [ux uy] = meshgrid(u(ii,:));
            r = dx*sum(Q.*fS(ux,uy),2); % (kron(Q,W).*FS)*1
            rhs(ii,:) = r;
        end
        
        % y coordinate
        for ii = 1:2*N+1
            [ux uy] = meshgrid(u(:,ii));
            r = dx*sum(Q.*fS(ux,uy),2);
            rhs(:,ii) = rhs(:,ii) + r;
        end                
        
        tau = .1*dx;
        rhs = -(rhs + tau*dx*(K*u + u*K'))/dx^2; % inv(kron(W,W))*rhs
        
        res = rk4a(INTRK)*res + dt*rhs;
        u   = u  + rk4b(INTRK)*res;        
    end
    
    if mod(i,interval)==0 || i==Nsteps
        Usnap(:,sk) = u(:);
        sk = sk + 1;
    end
    
    if mod(i,10)==0        
        surf(x,y,u)
        view(2)
        shading interp
        title(sprintf('step = %d / %d\n',i,Nsteps))
        drawnow
    end
end

