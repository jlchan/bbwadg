% SBP spectral method

clear
N = 50; % must be even
FinalTime = 2;

Np1D = N;

dx = 2/N;
h = 2*pi/N; % original h on [-pi,pi]
x = linspace(-1,1,N+1); x = x(1:end-1)';
column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)];
D = toeplitz(column,column([1 N:-1:2]))*pi;
W = dx*eye(Np1D);

Q = W*D;
K = D'*W*D;

%% make 2D

[x y] = meshgrid(x);
Np = Np1D^2;

%%

% Burgers
fS = @(uL,uR) (uL.^2 + uL.*uR + uR.^2)/6;
fD = @(uL,uR) max(abs(uL),abs(uR)).*(uR-uL);

%%


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
        for ii = 1:Np1D
            [ux uy] = meshgrid(u(ii,:));
            r = sum(Q.*fS(ux,uy),2); % (kron(Q,W).*FS)*1
            rhs(ii,:) = dx*r;
        end
        
        % y coordinate
        for ii = 1:Np1D
            [ux uy] = meshgrid(u(:,ii));
            r = sum(Q.*fS(ux,uy),2);
            rhs(:,ii) = rhs(:,ii) + dx*r;
        end                
        
        tau = .25*dx;
        rhs = -(rhs + tau*dx*(K*u + u*K))/dx^2; % inv(kron(W,W))*rhs
        
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

