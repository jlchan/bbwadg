% clear
K = 100;
%FinalTime = .35;
FinalTime = .25;

xv = linspace(-1,1,K+1)';
x = .5*(xv(1:end-1)+xv(2:end));
dx = 1/K;

CFL = .5;
tau = .1*dx;
Nmodes = 20;
snapshot_file = 'Usnap_euler2d_wall';

e = ones(K-1,1);
S = diag(e,1)-diag(e,-1);
S(1,:) = 0; S(:,1) = 0;
S(end,:) = 0; S(:,end) = 0;
S(1,2) = 1; S(2,1) = -1;
S(K-1,K) = 1; S(K,K-1) = -1;
% S(1,end) = -1; S(end,1) = 1;
S(1,1) = -1; S(end,end) = 1; % for BCs

S = sparse(S);
% KSx = sparse(2*eye(K) - abs(S))/dx;
KS = sparse(2*eye(K) - diag(e,1) - diag(e,-1))/dx;
KS(1,1) = KS(1,1)/2;
KS(K,K) = KS(K,K)/2;

% make 2D
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

% rho = ones(size(x)) + exp(-50*((x).^2+(y+.5).^2));
% % rho = 1 + (abs(x)<.5).*(abs(y)<.5);
% % xc = x + .5; yc = y + .5;
% % rho = 2*ones(size(x)) + .5*exp(-100*(xc.^2+0*yc.^2)); %.*sin(pi*xc).*sin(pi*yc);
% u = zeros(size(x));
% v = zeros(size(x));
% p = rho.^(gamma);
% 
% rhou = rho.*u;
% rhov = rho.*v;
% E    = p/(gamma-1) + .5*rho.*(u.^2+v.^2);

% % vortex
% p = ones(size(x));
% width = .1;
% sig = 1000;
% ff = 1-(1./(1+exp(-sig*(y-width))) + 1./(1+exp(sig*(y+width))));
% dv = .1*sin(2*pi*x).*ff;
% du = .5*ff;
% v = zeros(size(x)) + dv;
% u = -.5*ones(size(x)) + du;
% rho = ones(size(x));
% rhou = rho.*u;
% rhov = rho.*v;
% E    = p/(gamma-1) + .5*rho.*(u.^2+v.^2);

% % KH instab
% a = .1; sig = .1; %5*sqrt(2)*1e-3;
% ff = 1./(1+exp(-(y+.5)./(sig.^2)))-1./(1+exp(-(y-.5)./(sig.^2)));
% rho = ff + 1;
% v = a*sin(4*pi*x).*(exp(-(y-.5).^2/((sig/2).^2))+exp(-(y+.5).^2/((sig/2).^2)));
% u = ff - .5;
% p = 2.5*ones(size(x));
% rhou = rho.*u;
% rhov = rho.*v;
% E    = p/(gamma-1) + .5*rho.*(u.^2+v.^2);
% % surf(x,y,rho);return
%%

load(snapshot_file)
% load Usnap_euler2d_wall_novisc

U0 = reshape(Usnap(:,1),K^2,4);
u1 = Usnap(1:K^2,:);
u2 = Usnap((1:K^2) + K^2,:);
u3 = Usnap((1:K^2) + 2*K^2,:);
u4 = Usnap((1:K^2) + 3*K^2,:);

[Vr, Sr, ~] = svd([u1 u2 u3 u4 ...
    V1(u1,u2,u3,u4) V2(u1,u2,u3,u4) V3(u1,u2,u3,u4) V4(u1,u2,u3,u4)],0); 
sig = diag(Sr); 

% figure(1)
% semilogy(sig,'--','linewidth',2,'markersize',10)
% hold on
% grid on
% set(gca,'fontsize',16)

Vr = Vr(:,1:Nmodes);
Vrp = Vr;

Pr = Vr';

tol = sum(sig(Nmodes+1:end).^2)/sum(sig.^2)

rhoN = Vr'*U0(:,1);
rhouN = Vr'*U0(:,2);
rhovN = Vr'*U0(:,3);
EN = Vr'*U0(:,4);

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
         
     
dt = CFL*dx;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

res1 = zeros(size(rhoN));
res2 = zeros(size(rhoN));
res3 = zeros(size(rhoN));
res4 = zeros(size(rhoN));

rhs1 = zeros(size(rhoN));
rhs2 = zeros(size(rhoN));
rhs3 = zeros(size(rhoN));
rhs4 = zeros(size(rhoN));


sk = 2;
for i = 1:Nsteps
    for INTRK = 1:5                        
        
        rho  = Vr*rhoN;
        rhou = Vr*rhouN;
        rhov = Vr*rhovN;
        E    = Vr*EN;        
        
        v1 = Pr*V1(rho,rhou,rhov,E);
        v2 = Pr*V2(rho,rhou,rhov,E);
        v3 = Pr*V3(rho,rhou,rhov,E);
        v4 = Pr*V4(rho,rhou,rhov,E);        
        v1N = Vr*v1;
        v2N = Vr*v2;
        v3N = Vr*v3;
        v4N = Vr*v4;
        
        rho  = reshape(U1(v1N,v2N,v3N,v4N),K,K);
        rhou = reshape(U2(v1N,v2N,v3N,v4N),K,K);
        rhov = reshape(U3(v1N,v2N,v3N,v4N),K,K);
        E    = reshape(U4(v1N,v2N,v3N,v4N),K,K);        
        
        u = rhou./rho;
        v = rhov./rho;                        
            
        rhoL = [rho(:,1) rho];
        rhoR = [rho rho(:,end)];
        uL = [-u(:,1) u];
        uR = [u -u(:,end)];
        vL = [v(:,1) v];
        vR = [v v(:,end)];
        EL = [E(:,1) E];
        ER = [E E(:,end)];

        rhs1 = dx*diff(fxS1(rhoL,uL,vL,EL,rhoR,uR,vR,ER),[],2);
        rhs2 = dx*diff(fxS2(rhoL,uL,vL,EL,rhoR,uR,vR,ER),[],2);
        rhs3 = dx*diff(fxS3(rhoL,uL,vL,EL,rhoR,uR,vR,ER),[],2);
        rhs4 = dx*diff(fxS4(rhoL,uL,vL,EL,rhoR,uR,vR,ER),[],2);                        
        
        rhoL = [rho(1,:); rho];
        rhoR = [rho; rho(end,:)];
        uL = [u(1,:); u];
        uR = [u; u(end,:)];
        vL = [-v(1,:); v];
        vR = [v; -v(end,:)];
        EL = [E(1,:); E];
        ER = [E; E(end,:)];
        
        rhs1 = rhs1 + dx*diff(fyS1(rhoL,uL,vL,EL,rhoR,uR,vR,ER),[],1);
        rhs2 = rhs2 + dx*diff(fyS2(rhoL,uL,vL,EL,rhoR,uR,vR,ER),[],1);
        rhs3 = rhs3 + dx*diff(fyS3(rhoL,uL,vL,EL,rhoR,uR,vR,ER),[],1);
        rhs4 = rhs4 + dx*diff(fyS4(rhoL,uL,vL,EL,rhoR,uR,vR,ER),[],1);        
        
        % add dissipation                             
        rhs1 = rhs1 + tau*dx*(KS*rho + rho*KS);
        rhs2 = rhs2 + tau*dx*(KS*rhou + rhou*KS);
        rhs3 = rhs3 + tau*dx*(KS*rhov + rhov*KS);
        rhs4 = rhs4 + tau*dx*(KS*E + E*KS);
        
        % mass matrix = identity
        rhs1 = -Pr*rhs1(:)/dx^2;
        rhs2 = -Pr*rhs2(:)/dx^2;
        rhs3 = -Pr*rhs3(:)/dx^2;
        rhs4 = -Pr*rhs4(:)/dx^2;
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        res4 = rk4a(INTRK)*res4 + dt*rhs4;
        
        rhoN  = rhoN   + rk4b(INTRK)*res1;
        rhouN = rhouN  + rk4b(INTRK)*res2;
        rhovN = rhovN  + rk4b(INTRK)*res3;
        EN    = EN     + rk4b(INTRK)*res4;

    end       
    
    if mod(i,5) == 0
        clf
        surf(x,y,reshape(Vr*rhoN,K,K))
        shading interp
        view(2)
        
        Usnapi = reshape(Usnap(:,i),K*K,4);
        Uhi = Vrp*[rhoN rhouN rhovN EN];
        e = Uhi - Usnapi;
        err = sqrt(sum(sum(dx^2*(e.^2))))/sqrt(sum(sum(dx^2*Uhi.^2)));
        title(sprintf('time = %f, step %d / %d, err = %g\n',dt*i,i,Nsteps,err));
        
        drawnow
    end
end
