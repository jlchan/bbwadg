% SBP spectral method

clear

N = 40; % must be even
FinalTime = 5;

Np1D = N;

x = linspace(-1,1,N+1); x = x(1:end-1)';
% x = linspace(0,1,N+1); x = x(1:end-1)';

dx = (max(x)-min(x))/N;
h = 2*pi/N; % original h on [-pi,pi]

column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)];
D = toeplitz(column,column([1 N:-1:2]))*pi;
W = dx*eye(Np1D);

Q = W*D;
KS = D'*W*D;

FFT = zeros(N);
iFFT = zeros(N);
ei = zeros(N,1);
for i = 1:N
    ei(i) = 1;
    FFT(:,i) = fft(ei);
    iFFT(:,i) = ifft(ei);
    ei(i) = 0;
end

sfilter = ones(N,1);
kN = [0:N-1]/N;
Nc = .33;
alpha = 1;
p = 4;
sfilter(kN > Nc) = 0; %exp(-alpha*((kN(kN>Nc)-Nc)./(1-Nc)).^p);
F = real(iFFT*diag(sfilter)*FFT);

f = abs(x)<.5;
% f = tanh(10*sin(pi*(x-.5))); 
% f = exp(sin(pi*x));
f = exp(-500*x.^2);

plot(x,f,'o');
hold on
plot(x,D*f,'s-');
plot(x,F*D*f,'*--');
% plot(x,D*f,'o--');


%% make 2D

[x y] = meshgrid(x);
Np = Np1D^2;

%% Euler 2D

gamma = 1.4;

rhoe = @(rho,rhou,rhov,E) E - .5*(rhou.^2+rhov.^2)./rho;
pcons = @(rho,rhou,rhov,E) (gamma-1)*rhoe(rho,rhou,rhov,E);
sfun = @(rho,rhou,rhov,E) log((gamma-1)*rhoe(rho,rhou,rhov,E)./(rho.^gamma));
V1 = @(rho,rhou,rhov,E) (-E + rhoe(rho,rhou,rhov,E).*(gamma + 1 - sfun(rho,rhou,rhov,E)))./(rhoe(rho,rhou,rhov,E));
V2 = @(rho,rhou,rhov,E) rhou./(rhoe(rho,rhou,rhov,E));
V3 = @(rho,rhou,rhov,E) rhov./(rhoe(rho,rhou,rhov,E));
V4 = @(rho,rhou,rhov,E) (-rho)./(rhoe(rho,rhou,rhov,E));

sV = @(V1,V2,V3,V4) gamma - V1 + (V2.^2+V3.^2)./(2*V4);
rhoeV  = @(V1,V2,V3,V4) ((gamma-1)./((-V4).^gamma)).^(1/(gamma-1)).*exp(-sV(V1,V2,V3,V4)/(gamma-1));
U1 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(-V4);
U2 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(V2);
U3 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(V3);
U4 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(1-(V2.^2+V3.^2)./(2*V4));

avg = @(x,y) .5*(x+y);
pfun = @(rho,u,v,E) (gamma-1)*(E-.5*rho.*(u.^2+v.^2));
beta = @(rho,u,v,E) rho./(2*pfun(rho,u,v,E));
pavg     = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER)     avg(rhoL,rhoR)./(2*avg(beta(rhoL,uL,vL,EL),beta(rhoR,uR,vR,ER)));
plogmean = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR)./(2*logmean(beta(rhoL,uL,vL,EL),beta(rhoR,uR,vR,ER)));
vnormavg = @(uL,vL,uR,vR) 2*(avg(uL,uR).^2 + avg(vL,vR).^2) - (avg(uL.^2,uR.^2) + avg(vL.^2,vR.^2));

fxS1 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(uL,uR);
fxS2 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(uL,uR).^2 + pavg(rhoL,uL,vL,EL,rhoR,uR,vR,ER);
fxS3 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(uL,uR).*avg(vL,vR);
fxS4 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) (plogmean(rhoL,uL,vL,EL,rhoR,uR,vR,ER)/(gamma-1) + pavg(rhoL,uL,vL,EL,rhoR,uR,vR,ER) + .5*logmean(rhoL,rhoR).*vnormavg(uL,vL,uR,vR)).*avg(uL,uR);

fyS1 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(vL,vR);
fyS2 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(uL,uR).*avg(vL,vR);
fyS3 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(vL,vR).^2 + pavg(rhoL,uL,vL,EL,rhoR,uR,vR,ER);
fyS4 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) (plogmean(rhoL,uL,vL,EL,rhoR,uR,vR,ER)/(gamma-1) + pavg(rhoL,uL,vL,EL,rhoR,uR,vR,ER) + .5*logmean(rhoL,rhoR).*vnormavg(uL,vL,uR,vR)).*avg(vL,vR);

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

%%

% Munz IC
p = ones(size(x));
width = .1;
sig = 1000;
ff = 1-(1./(1+exp(-sig*(y-width))) + 1./(1+exp(sig*(y+width))));
dv = .1*sin(2*pi*x).*ff;
du = .5*ff;
v = zeros(size(x)) + dv;
u = -.5*ones(size(x)) + du;

rho = ones(size(x));
rhou = rho.*u;
rhov = rho.*v;
E    = p/(gamma-1) + .5*rho.*(u.^2+v.^2);

% % pulse
% a = 150;
% r2 = (x-.5).^2 + (y-.5).^2;
% rho = 1+exp(-a*r2);
% rhou = 0*rho;
% rhov = 0*rho;
% E = rho.^(gamma);

% KH instab
a = .1; sig = .1; %5*sqrt(2)*1e-3;
ff = 1./(1+exp(-(y-.25)./(sig.^2)))-1./(1+exp(-(y-.75)./(sig.^2)));
rho = ff + 1;
% sig = 5*sqrt(2)*1e-3;
v = a*sin(4*pi*x).*(exp(-(y-.25).^2/((sig/2).^2))+exp(-(y-.75).^2/((sig/2).^2)));
u = ff - .5;
p = 2.5*ones(size(x));

rhou = rho.*u;
rhov = rho.*v;
E    = p/(gamma-1) + .5*rho.*(u.^2+v.^2);

% surf(x,y,v); xlabel('x'); return

%%

dt = .5*dx;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

interval = inf;
Usnap = zeros(4*Np,ceil(Nsteps/interval)+1);
Usnap(:,1) = [rho(:);rhou(:);rhov(:);E(:)];

res1 = zeros(size(u));
res2 = zeros(size(u));
res3 = zeros(size(u));
res4 = zeros(size(u));

rhs1 = zeros(size(u));
rhs2 = zeros(size(u));
rhs3 = zeros(size(u));
rhs4 = zeros(size(u));

sk = 2;
for i = 1:Nsteps
    for INTRK = 1:5
        
        u = rhou./rho;
        v = rhov./rho;
        b = beta(rho,u,v,E);
        
        % x coordinate
        for ii = 1:Np1D
            [rhox rhoy] = meshgrid(rho(ii,:));
            [ux uy] = meshgrid(u(ii,:));
            [vx vy] = meshgrid(v(ii,:));
            [betax betay] = meshgrid(b(ii,:));
            
            rholog = logmean(rhox,rhoy);
            rhoavg = avg(rhox,rhoy);
            uavg = avg(ux,uy);
            vavg = avg(vx,vy);
            vnavg = 2*(uavg.^2 + vavg.^2) - (avg(ux.^2,uy.^2) + avg(vx.^2,vy.^2));
            pa = rhoavg./(2*avg(betax,betay));
            
            FxS1 = rholog.*uavg;      
            FxS2 = FxS1.*uavg + pa;   
            FxS3 = FxS1.*vavg;              
            f4aux = rholog./(2*(gamma-1)*logmean(betax,betay)) + pa + .5*rholog.*vnavg;
            FxS4 = f4aux.*uavg;                        

            rhs1(ii,:) = dx*sum(Q.*FxS1,2);
            rhs2(ii,:) = dx*sum(Q.*FxS2,2);
            rhs3(ii,:) = dx*sum(Q.*FxS3,2);
            rhs4(ii,:) = dx*sum(Q.*FxS4,2);
        end
        
        % y coordinate
        for ii = 1:Np1D
            [rhox rhoy] = meshgrid(rho(:,ii));
            [ux uy] = meshgrid(u(:,ii));
            [vx vy] = meshgrid(v(:,ii));
            [betax betay] = meshgrid(b(:,ii));
            
            rholog = logmean(rhox,rhoy);
            rhoavg = avg(rhox,rhoy);
            uavg = avg(ux,uy);
            vavg = avg(vx,vy);
            vnavg = 2*(uavg.^2 + vavg.^2) - (avg(ux.^2,uy.^2) + avg(vx.^2,vy.^2));
            pa = rhoavg./(2*avg(betax,betay));
            
            FyS1 = rholog.*vavg;
            FyS2 = FyS1.*uavg;
            FyS3 = FyS1.*vavg + pa;
            f4aux = rholog./(2*(gamma-1)*logmean(betax,betay)) + pa + .5*rholog.*vnavg;
            FyS4 = f4aux.*vavg;
                        
            rhs1(:,ii) = rhs1(:,ii) + dx*sum(Q.*FyS1,2);
            rhs2(:,ii) = rhs2(:,ii) + dx*sum(Q.*FyS2,2);
            rhs3(:,ii) = rhs3(:,ii) + dx*sum(Q.*FyS3,2);
            rhs4(:,ii) = rhs4(:,ii) + dx*sum(Q.*FyS4,2);
        end
        
        tau = 1e-4;
        rhs1 = -(rhs1 + tau*dx*(KS*rho + rho*KS))/dx^2;
        rhs2 = -(rhs2 + tau*dx*(KS*rhou + rhou*KS))/dx^2;
        rhs3 = -(rhs3 + tau*dx*(KS*rhov + rhov*KS))/dx^2;
        rhs4 = -(rhs4 + tau*dx*(KS*E + E*KS))/dx^2;
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        res4 = rk4a(INTRK)*res4 + dt*rhs4;
        
        rho  = rho  + rk4b(INTRK)*res1;
        rhou = rhou + rk4b(INTRK)*res2;
        rhov = rhov + rk4b(INTRK)*res3;
        E    = E    + rk4b(INTRK)*res4;        
    end
    
    if mod(i,interval)==0 || i==Nsteps
        Usnap(:,sk) = [rho(:);rhou(:);rhov(:);E(:)];
        sk = sk + 1;
    end
    
    if (mod(i,10)==0 || i==Nsteps)
        surf(x,y,rho)
        view(2)
        colorbar
        xlabel('x')
        ylabel('y')
%         axis([-1,1,-1,1,-1,1])
        shading interp
        title(sprintf('step = %d / %d, time = %g\n',i,Nsteps,dt*i))
        drawnow
    end
end

