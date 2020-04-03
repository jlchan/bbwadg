clear
K = 1000;
FinalTime = .7;
% FinalTime = .2;
xv = linspace(0,1,K+1)';
x = .5*(xv(1:end-1)+xv(2:end));
dx = 1/K;

tau = .5*dx;       
% tau = 5e-4; 
% tau = .25;

e = ones(K-1,1);
S = diag(e,1)-diag(e,-1);
S(1,:) = 0; S(:,1) = 0;
S(end,:) = 0; S(:,end) = 0;
S(1,end) = -1; S(end,1) = 1;
S(1,2) = 1; S(2,1) = -1;
S(K-1,K) = 1; S(K,K-1) = -1;

S = sparse(S);
%KS = sparse(2*eye(K) - abs(S))/dx^2;

e = ones(K,1);
D = diag(e) - diag(e(2:end),1);
D = D(1:end-1,:);
KS = sparse(D'*D)/dx;

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

c = @(rho,u,E) sqrt(gamma*pfun(rho,u,E)./rho);
lam = @(rho,u,E) abs(u) + c(rho,u,E);
llf = @(rhoL,rhoR,uL,uR,EL,ER) .5*(lam(rhoL,uL,EL)+lam(rhoR,uR,ER));
fD1 = @(rhoL,rhoR,uL,uR,EL,ER) llf(rhoL,rhoR,uL,uR,EL,ER).*(rhoL-rhoR);
fD2 = @(rhoL,rhoR,uL,uR,EL,ER) llf(rhoL,rhoR,uL,uR,EL,ER).*(rhoL.*uL-rhoR.*uR);
fD3 = @(rhoL,rhoR,uL,uR,EL,ER) llf(rhoL,rhoR,uL,uR,EL,ER).*(EL-ER);

%% set initial cond

opt = 2;

if opt==0    
    % sine solution
    t = 0;    
    rhoex = @(x) (2 + sin(pi*((2*x-1) - t)));
    
    a = .25;
    uex = @(x) ones(size(x));
    pex = @(x) ones(size(x)) + a*exp(-25*(2*x-1).^2);    
    
elseif opt==1
    % afkam, hesthaven    
    rhoex = @(x) .5 + .2*cos(2*pi*x);
    uex = @(x) 1.5 + 0*x;
    pex = @(x) .5 + .2*sin(2*pi*x);
elseif opt==2
    
    rhoex = @(x) 2 + exp(-25*(2*x-1).^2);
    uex = @(x) 0*x;
    pex = @(x) 1+0*x;
    pex = @(x) rhoex(x).^gamma;
    
end

mex = @(x) rhoex(x).*uex(x);
Eex = @(x) pex(x)./(gamma-1) + .5*rhoex(x).*uex(x).^2;

rho = rhoex(x);
m = mex(x);
E = Eex(x);

%%
trainROM = 0;
if trainROM
    load Usnap_euler
    
    U0 = reshape(Usnap(:,1),K,3);
    
    Nsample = 3; %ceil(size(Usnap,2)/500); % downsample
    Us1 = Usnap((1:K),1:Nsample:end);
    Us2 = Usnap((1:K)+K,1:Nsample:end);
    Us3 = Usnap((1:K)+2*K,1:Nsample:end);
    
    % add entropy variables to snapshots
    Us = [Us1 Us2 Us3 V1(Us1,Us2,Us3) V2(Us1,Us2,Us3) V3(Us1,Us2,Us3)];
    [Vr,Sr,~] = svd(Us,0);
    
    Nmodes = 35;
    
    Vr = Vr(:,1:Nmodes);    
    Mr = Vr'*Vr*dx;    
    Pr = Mr\(Vr'*dx);    
    
else
    Vr = speye(K);
    Pr = speye(K);
end

rho = Vr'*rho;
m = Vr'*m;
E = Vr'*E;

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
res1 = zeros(size(rho));
res2 = zeros(size(rho));
res3 = zeros(size(rho));


interval = 1;
if trainROM
    Fsnap = zeros(3*K,ceil(Nsteps/interval));
    sk = 1;
else
    Usnap = zeros(3*K,ceil(Nsteps/interval));
    Usnap(:,1) = [rho;m;E];
    sk = 2;
end
for i = 1:Nsteps
    for INTRK = 1:5   
        if trainROM
            rhoq = Vr*rho;
            mq = Vr*m;
            Eq = Vr*E;
            v1 = Vr*Pr*V1(rhoq,mq,Eq);
            v2 = Vr*Pr*V2(rhoq,mq,Eq);
            v3 = Vr*Pr*V3(rhoq,mq,Eq);
            rhoq = U1(v1,v2,v3);
            mq = U2(v1,v2,v3);
            Eq = U3(v1,v2,v3);
        else
            rhoq = rho;
            mq = m;
            Eq = E;
        end
        
        uq = mq./rhoq;        
                
        rhoL = [rhoq(end) rhoq(:)'];
        rhoR = [rhoq(:)' rhoq(1)];        
        uL = [uq(end) uq(:)'];
        uR = [uq(:)' uq(1)];        
        EL = [Eq(end) Eq(:)'];
        ER = [Eq(:)' Eq(1)];
                
        r1 = diff(fS1(rhoL,rhoR,uL,uR,EL,ER))';
        r2 = diff(fS2(rhoL,rhoR,uL,uR,EL,ER))';
        r3 = diff(fS3(rhoL,rhoR,uL,uR,EL,ER))';
                
        if INTRK==5
           rhstest(i) = full(sum(sum(V1(rhoq,mq,Eq).*r1+V2(rhoq,mq,Eq).*r2+V3(rhoq,mq,Eq).*r3)));
        end                
        
        % add dissipation + invert mass
        d1 = KS*rhoq;
        d2 = KS*mq;
        d3 = KS*Eq;
        
%         d1 = diff(fD1(rhoL,rhoR,uL,uR,EL,ER))'/dx;
%         d2 = diff(fD2(rhoL,rhoR,uL,uR,EL,ER))'/dx;
%         d3 = diff(fD3(rhoL,rhoR,uL,uR,EL,ER))'/dx;
        
        rhs1 = -Pr*((r1 + tau*d1)/dx);
        rhs2 = -Pr*((r2 + tau*d2)/dx);
        rhs3 = -Pr*((r3 + tau*d3)/dx);
        
        if INTRK==5
           dtest(i) = full(sum(sum(V1(rhoq,mq,Eq).*d1+V2(rhoq,mq,Eq).*d2+V3(rhoq,mq,Eq).*d3)));
        end                
%         tau = .25;
%         d1 = diff(fD1(rhoL,rhoR,uL,uR,EL,ER))';
%         d2 = diff(fD2(rhoL,rhoR,uL,uR,EL,ER))';
%         d3 = diff(fD3(rhoL,rhoR,uL,uR,EL,ER))';
%         rhs1 = -Pr*((r1 + tau*d1)/dx);
%         rhs2 = -Pr*((r2 + tau*d2)/dx);
%         rhs3 = -Pr*((r3 + tau*d3)/dx);

        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        rho = rho + rk4b(INTRK)*res1;
        m   = m   + rk4b(INTRK)*res2;
        E   = E   + rk4b(INTRK)*res3;
    end      
    
    
    Sq = -rho.*s(rho,m,E);
    fom_entropy(i) = sum(sum(dx.*Sq));
    
    if trainROM
        if mod(i,interval)==0
            Fsnap(:,sk) = [r1;r2;r3];
            sk = sk + 1;
        end
    else
        if mod(i,interval)==0
            Usnap(:,sk) = [rho;m;E];
            sk = sk + 1;
        end
    end
    if mod(i,10) == 0
%         plot(x,d1,'o')
        plot(x,Vr*rho,'o')
%         plot(x,pfun(rho,m./rho,E),'o')
%         hold on
%         plot(x,m./rho,'o')
%         plot(x,E,'o')
%         hold off
%         axis([0,1,.8,4])  
        %plot(x,pfun(rho,m./rho,E),'o--')
%         axis([-1,1,.8,1.25])  
        title(sprintf('time = %f, step %d / %d, rhstest = %g, dtest = %g\n',dt*i,i,Nsteps,rhstest(i),dtest(i)));
        drawnow
    end
end

figure
semilogy(dt*(1:Nsteps),fom_entropy,'--','linewidth',2)

