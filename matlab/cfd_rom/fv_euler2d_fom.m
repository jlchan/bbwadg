clear
K = 200;
FinalTime = .5;
xv = linspace(-1,1,K+1)';
x = .5*(xv(1:end-1)+xv(2:end));
dx = 1/K;

e = ones(K-1,1);
S = diag(e,1)-diag(e,-1);
S(1,:) = 0; S(:,1) = 0;
S(end,:) = 0; S(:,end) = 0;
S(1,end) = -1; S(end,1) = 1;
S(1,2) = 1; S(2,1) = -1;
S(K-1,K) = 1; S(K,K-1) = -1;

S = sparse(S);
KSx = sparse(2*eye(K) - abs(S))/dx;
KSy = sparse(2*eye(K) - diag(e,1) - diag(e,-1))/dx;
KSy(1,1) = KSy(1,1)/2;
KSy(K,K) = KSy(K,K)/2;

%% make 2D

[x y] = meshgrid(x);


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

%% set initial cond

rho = ones(size(x)) + exp(-100*(x.^2+(y+.5).^2));
p = rho.^(gamma);
u = zeros(size(x));
v = zeros(size(x));
rhou = rho.*u;
rhov = rho.*v;
E    = p/(gamma-1) + .5*rho.*(u.^2+v.^2);

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
         
     
dt = .75*dx;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

res1 = zeros(size(x));
res2 = zeros(size(x));
res3 = zeros(size(x));
res4 = zeros(size(x));

rhs1 = zeros(size(x));
rhs2 = zeros(size(x));
rhs3 = zeros(size(x));
rhs4 = zeros(size(x));

interval = inf;
Usnap = zeros(4*K*K,ceil(Nsteps/interval));
Usnap(:,1) = [rho(:);rhou(:);rhov(:);E(:)];
sk = 2;
for i = 1:Nsteps
    for INTRK = 1:5                        
        
        u = rhou./rho;
        v = rhov./rho;
                
        tau = .1*dx;
            
        rhoL = [rho(:,end) rho];
        rhoR = [rho rho(:,1)];
        uL = [u(:,end) u];
        uR = [u u(:,1)];
        vL = [v(:,end) v];
        vR = [v v(:,1)];
        EL = [E(:,end) E];
        ER = [E E(:,1)];
        
        rhs1 = diff(fxS1(rhoL,uL,vL,EL,rhoR,uR,vR,ER),[],2);
        rhs2 = diff(fxS2(rhoL,uL,vL,EL,rhoR,uR,vR,ER),[],2);
        rhs3 = diff(fxS3(rhoL,uL,vL,EL,rhoR,uR,vR,ER),[],2);
        rhs4 = diff(fxS4(rhoL,uL,vL,EL,rhoR,uR,vR,ER),[],2);                        
        
        rhoL = [rho(1,:); rho];
        rhoR = [rho; rho(end,:)];
        uL = [u(1,:); u];
        uR = [u; u(end,:)];
        vL = [-v(1,:); v];
        vR = [v; -v(end,:)];
        EL = [E(1,:); E];
        ER = [E; E(end,:)];
        
        rhs1 = rhs1 + diff(fyS1(rhoL,uL,vL,EL,rhoR,uR,vR,ER),[],1);
        rhs2 = rhs2 + diff(fyS2(rhoL,uL,vL,EL,rhoR,uR,vR,ER),[],1);
        rhs3 = rhs3 + diff(fyS3(rhoL,uL,vL,EL,rhoR,uR,vR,ER),[],1);
        rhs4 = rhs4 + diff(fyS4(rhoL,uL,vL,EL,rhoR,uR,vR,ER),[],1);        
        
        % add dissipation
        rhs1 = rhs1 + tau*(KSy*rho + rho*KSx);
        rhs2 = rhs2 + tau*(KSy*rhou + rhou*KSx);
        rhs3 = rhs3 + tau*(KSy*rhov + rhov*KSx);
        rhs4 = rhs4 + tau*(KSy*E + E*KSx);
        
        rhs1 = -rhs1/dx;
        rhs2 = -rhs2/dx;
        rhs3 = -rhs3/dx;
        rhs4 = -rhs4/dx;
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        res4 = rk4a(INTRK)*res4 + dt*rhs4;
        
        rho  = rho  + rk4b(INTRK)*res1;
        rhou = rhou  + rk4b(INTRK)*res2;
        rhov = rhov  + rk4b(INTRK)*res3;
        E    = E  + rk4b(INTRK)*res4;

    end
    
    if mod(i,interval)==0
        Usnap(:,sk) = [rho(:);rhou(:);rhov(:);E(:)];
        sk = sk + 1;
    end
    
    if mod(i,1) == 0
        surf(x,y,reshape(rho,K,K))
        shading interp
        title(sprintf('time = %f, step %d / %d\n',dt*i,i,Nsteps));
        drawnow
    end
end