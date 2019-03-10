clear
K = 100;
FinalTime = 4;
xv = linspace(-1,1,K+1)';
x = .5*(xv(1:end-1)+xv(2:end));
dx = 2/K;

e = ones(K-1,1);
S = diag(e,1)-diag(e,-1);
S(1,:) = 0; S(:,1) = 0;
S(end,:) = 0; S(:,end) = 0;
S(1,end) = -1; S(end,1) = 1;
S(1,2) = 1; S(2,1) = -1;
S(K-1,K) = 1; S(K,K-1) = -1;

S = sparse(S);
KS = sparse(2*eye(K) - abs(S))/dx;
SD = abs(S); % defunct?

%% Burgers

fS = @(uL,uR) (uL.^2 + uL.*uR + uR.^2)/6;
fD = @(uL,uR) max(abs(uL),abs(uR)).*(uR-uL);
fD = @(uL,uR) (uR-uL);
% fD = @(uL,uR) 0*max(abs(uL),abs(uR)).*(uR-uL);

% fS = @(uL,uR) (uL + uR)/2;
% fD = @(uL,uR) (uR-uL);

%% Euler fluxes

gamma = 1.4;
rhoe = @(rho,m,E) E - .5*m.^2./rho;
s = @(rho,m,E) log((gamma-1).*rhoe(rho,m,E)./(rho.^gamma));

V1 = @(rho,m,E) (-E + rhoe(rho,m,E).*(gamma + 1 - s(rho,m,E)))./(rhoe(rho,m,E));
V2 = @(rho,m,E) (m)./(rhoe(rho,m,E));
V3 = @(rho,m,E) (-rho)./(rhoe(rho,m,E));

sV = @(V1,V2,V3) gamma - V1 + V2.^2./(2*V3);
rhoeV  = @(V1,V2,V3) ((gamma-1)./((-V3).^gamma)).^(1/(gamma-1)).*exp(-sV(V1,V2,V3)/(gamma-1));
U1 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(-V3);
U2 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(V2);
U3 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(1-V2.^2./(2*V3));

avg = @(uL,uR) .5*(uL+uR);
pfun = @(rho,u,E) (gamma-1)*(E-.5*rho.*u.^2);
beta = @(rho,u,E) rho./(2*pfun(rho,u,E));

pavg = @(rhoL,rhoR,uL,uR,EL,ER) avg(rhoL,rhoR)./(2*avg(beta(rhoL,uL,EL),beta(rhoR,uR,ER)));
fS1 = @(rhoL,rhoR,uL,uR,EL,ER) logmean(rhoL,rhoR).*avg(uL,uR);
fS2 = @(rhoL,rhoR,uL,uR,EL,ER) pavg(rhoL,rhoR,uL,uR,EL,ER) + avg(uL,uR).*fS1(rhoL,rhoR,uL,uR,EL,ER);
fS3 = @(rhoL,rhoR,uL,uR,EL,ER) fS1(rhoL,rhoR,uL,uR,EL,ER)...
    .*(1./(2*(gamma-1).*logmean(beta(rhoL,uL,EL),beta(rhoR,uR,ER))) - .5*avg(uL.^2,uR.^2)) ...
    + avg(uL,uR).*fS2(rhoL,rhoR,uL,uR,EL,ER);


%% set initial cond

opt = 0;

if opt==0    
    % sine solution
    t = 0;
    rhoex = @(x) (2 + sin(pi*(x - t)));
    
    a = .25;
    uex = @(x) ones(size(x));
    pex = @(x) ones(size(x)) + a*exp(-25*x.^2);    
        
elseif opt==1
    
    % pulse condition    
    rhoex = @(x) 1+exp(-100*x.^2);
    pex = @(x) rhoex(x).^gamma;
    uex = @(x) 0*x;
    
elseif opt==2
            
    % sod initial condition
    rhoL = 1; rhoR = .125;
    pL = 1; pR = .1;
    rhoex = @(x) (rhoL*(x < 0) + rhoR*(x > 0));
    pex = @(x) (pL*(x < 0) + pR*(x > 0));
    uex = @(x) 0*x;    
    
end

mex = @(x) rhoex(x).*uex(x);
Eex = @(x) pex(x)./(gamma-1) + .5*rhoex(x).*uex(x).^2;

rho = rhoex(x);
m = mex(x);
E = Eex(x);

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
         
     
dt = .5*dx;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;
res1 = zeros(size(x));
res2 = zeros(size(x));
res3 = zeros(size(x));

interval = 1;
Usnap = zeros(3*K,ceil(Nsteps/interval));
Usnap(:,1) = [rho;m;E];
sk = 2;
for i = 1:Nsteps
    for INTRK = 1:5                        
        u = m./rho;        
                
        [rhox rhoy] = meshgrid(rho);
        [ux uy] = meshgrid(u);
        [Ex Ey] = meshgrid(E);
        
        rhs1 = sum(S.*fS1(rhox,rhoy,ux,uy,Ex,Ey),2);
        rhs2 = sum(S.*fS2(rhox,rhoy,ux,uy,Ex,Ey),2);
        rhs3 = sum(S.*fS3(rhox,rhoy,ux,uy,Ex,Ey),2);
        
%         tau = .5;
%         cvel = sqrt(gamma*pfun(rho,u,E)./rho);
%         lm   = (abs(u) + cvel);
%         Lfc  = .5*(SD*lm); % average wavespeed
%         rhs1 = sum(S.*fS1(rhox,rhoy,ux,uy,Ex,Ey),2) + tau*KS*(Lfc.*rho);
%         rhs2 = sum(S.*fS2(rhox,rhoy,ux,uy,Ex,Ey),2) + tau*KS*(Lfc.*m);
%         rhs3 = sum(S.*fS3(rhox,rhoy,ux,uy,Ex,Ey),2) + tau*KS*(Lfc.*E);
        
        if INTRK==5
           rhstest(i) = full(sum(sum(V1(rho,m,E).*rhs1+V2(rho,m,E).*rhs2+V3(rho,m,E).*rhs3)));
        end                
        
        % add dissipation + invert mass
        tau = .1*dx;
        rhs1 = -(rhs1 + tau*KS*rho)/dx;
        rhs2 = -(rhs2 + tau*KS*m)/dx;
        rhs3 = -(rhs3 + tau*KS*E)/dx;

        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        rho = rho  + rk4b(INTRK)*res1;
        m   = m  + rk4b(INTRK)*res2;
        E   = E  + rk4b(INTRK)*res3;
    end      
    
    if mod(i,interval)==0
        Usnap(:,sk) = [rho;m;E];
        sk = sk + 1;
    end
    
    if mod(i,10) == 0
        plot(x,rho,'o--')
        axis([-1,1,.8,4])  
        %plot(x,pfun(rho,m./rho,E),'o--')
        %axis([-1,1,.8,1.25])  
        title(sprintf('rhstest = %g\n',rhstest(i)));
        drawnow
    end
end
