% SBP spectral method

clear

N = 40; % must be even
FinalTime = 5.0;

CFL = .25;

Lx = 20;
x = linspace(0,Lx,N+1); x = x(1:end-1)';
y = x; %linspace(0,10,N+1); y = y(1:end-1)';

dx = x(2)-x(1);
h = 2*pi/N; % original h on [-pi,pi]


column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]*2*pi/Lx;
D = toeplitz(column,column([1 N:-1:2]));
W = dx*eye(N);

Q = W*D;

% k = 1*2/Lx;
% plot(D*sin(k*pi*x)-k*pi*cos(k*pi*x),'o')
% hold on
% plot()

%% make 2D

[x y] = meshgrid(x,y);

%% Euler 2D

global gamma
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

% pulse
a = 1/2;
x0 = 10;
y0 = 10;
r2 = (x-x0).^2 + (y-y0).^2;
rho = 2+exp(-a*r2);
u = .5*ones(size(x));
v = 0*rho;
p = rho.^(gamma);

[rho u v p] = vortexSolution(x,y,0);

rhou = rho.*u;
rhov = rho.*v;
E    = p/(gamma-1) + .5*rho.*(u.^2+v.^2);

% surf(x,y,rhou)
% return

Nq = N;
[Vp,~,xp1D] = SincBasis(N,Nq);
xp = linspace(0,Lx,Nq+1);
xp = xp(1:end-1);
[xp yp] = meshgrid(xp);
dxp = xp(2)-xp(1);

%%


dt = CFL*dx;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

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
        for ii = 1:N
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
        for ii = 1:N
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
        

        
        if (INTRK==5)
            rhstest(i) = sum(sum(V1(rho,rhou,rhov,E).*rhs1 + ...
                V2(rho,rhou,rhov,E).*rhs2 + V3(rho,rhou,rhov,E).*rhs3 + V4(rho,rhou,rhov,E).*rhs4));
        end
        
        rhs1 = -(2*rhs1)/dx^2;
        rhs2 = -(2*rhs2)/dx^2;
        rhs3 = -(2*rhs3)/dx^2;
        rhs4 = -(2*rhs4)/dx^2;
       
        
%         tau = 1e-4;
%         rhs1 = -(rhs1 + tau*dx*(KS*rho + rho*KS))/dx^2;
%         rhs2 = -(rhs2 + tau*dx*(KS*rhou + rhou*KS))/dx^2;
%         rhs3 = -(rhs3 + tau*dx*(KS*rhov + rhov*KS))/dx^2;
%         rhs4 = -(rhs4 + tau*dx*(KS*E + E*KS))/dx^2;
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        res4 = rk4a(INTRK)*res4 + dt*rhs4;
                
        rho  = rho  + rk4b(INTRK)*res1;
        rhou = rhou + rk4b(INTRK)*res2;
        rhov = rhov + rk4b(INTRK)*res3;
        E    = E    + rk4b(INTRK)*res4;        
        
%         [rhs1(:) rhs2(:) rhs3(:) rhs4(:)]        
%         [res1(:) res2(:) res3(:) res4(:)]
%         [rho(:) rhou(:) rhov(:) E(:)]
%         keyboard

    end
        
    if (mod(i,5)==0 || i==Nsteps)
        surf(xp,yp,(Vp*rho)*Vp')
        view(2)
        colorbar
        xlabel('x')
        ylabel('y')
%         axis([-1,1,-1,1,-1,1])
        shading interp
        title(sprintf('step = %d / %d, time = %g, rhstest = %g\n',i,Nsteps,dt*i,rhstest(i)))
        drawnow
    end
end

[rhoex u v p] = vortexSolution(xp,yp,FinalTime);
rhouex = rhoex.*u;
rhovex = rhoex.*v;
Eex    = p/(gamma-1) + .5*rhoex.*(u.^2+v.^2);

err = sqrt(sum(sum(dx*dx*(((Vp*rho)*Vp'-rhoex).^2 + ((Vp*rhou)*Vp'-rhouex).^2 + ((Vp*rhov)*Vp'-rhovex).^2 + ((Vp*E)*Vp'-Eex).^2))))


subplot(2,1,1)
surf(xp,yp,(Vp*rho)*Vp')
view(2)
colorbar
xlabel('x')
ylabel('y')
shading interp
title(sprintf('Err = %g\n',err),'fontsize',14)
axis equal
subplot(2,1,2)
surf(xp,yp,rhoex-(Vp*rho)*Vp')
view(2)
colorbar
xlabel('x')
ylabel('y')
shading interp
axis equal

function [rho u v p] = vortexSolution(x,y,t)

global gamma
x0 = 5;
y0 = 10;
beta = 5;
r2 = (x-x0-t).^2 + (y-y0).^2;

u = 1 - beta*exp(1-r2).*(y-y0)/(2*pi);
v = beta*exp(1-r2).*(x-x0-t)/(2*pi);
rho = 1 - (1/(8*gamma*pi^2))*(gamma-1)/2*(beta*exp(1-r2)).^2;
rho = rho.^(1/(gamma-1));
p = rho.^gamma;

% sine sol
%k = [1 .5];
%rho = (2 + sin(pi*(1*(x-k(1)*t)+2*(y-k(2)*t))*2/20));
%u = ones(size(x))*k(1);
%v = ones(size(x))*k(2);
rho = 2 + sin(pi*(x-t)*2/20);
u = ones(size(x));
v = zeros(size(x));
p = ones(size(x));

% if 1
%     % pulse condition
%     x0 = 7.5;
%     rho = 2 + (abs(x-x0) < 2.5);
%     u = 0*rho;
%     v = 0*rho;
%     p = rho.^gamma;
%     
% end

end