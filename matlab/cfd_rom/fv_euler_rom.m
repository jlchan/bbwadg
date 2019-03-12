clear
K = 1000;
FinalTime = .7;
xv = linspace(0,1,K+1)';
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
KS = sparse(2*eye(K) - abs(S))/dx;

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

opt = 1;

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
    uex = @(x) 1.5;
    pex = @(x) .5 + .2*sin(2*pi*x);
end

mex = @(x) rhoex(x).*uex(x);
Eex = @(x) pex(x)./(gamma-1) + .5*rhoex(x).*uex(x).^2;

rho = rhoex(x);
m = mex(x);
E = Eex(x);

%% construct ROM

if 1
    load Usnap_euler
    
    U0 = reshape(Usnap(:,1),K,3);
    %  rho = Vr'*rho;
    %  m = Vr'*m;
    %  E = Vr'*E;    

    Nsample = 5; %ceil(size(Usnap,2)/500); % downsample
    Us1 = Usnap((1:K),1:Nsample:end);
    Us2 = Usnap((1:K)+K,1:Nsample:end);
    Us3 = Usnap((1:K)+2*K,1:Nsample:end);
    
    % add entropy variables to snapshots
    Us = [Us1 Us2 Us3 V1(Us1,Us2,Us3) V2(Us1,Us2,Us3) V3(Us1,Us2,Us3)];
%     Us = [Us1 Us2 Us3];
    [Vr,Sr,~] = svd(Us,0);
    
    Nmodes = 50;
    Vr = Vr(:,1:Nmodes);
    Vrp = Vr;
    
    sig = diag(Sr);
    tol = sqrt(sum(sig(Nmodes+1:end).^2)/sum(sig.^2))
    
    % add range + constants to snapshots
    [Vtest Stest, ~] = svd([Vr S*Vr]);
    sigt = diag(Stest);        
    Vtest = Vtest(:,sigt > 1e-8);
    Vtest = orth([ones(size(x)) Vtest]); % ensure 1
    
    Sfull = S;
    Kfull = KS;
        
    rho = Vr'*U0(:,1);
    m = Vr'*U0(:,2);
    E = Vr'*U0(:,3);        
    
    % empirical cubature
    if 1        
%         % target space
%         sk = 1;
%         for i = 1:size(Vtest,2)
%             for j = i:size(Vtest,2)
%                 Vmass(:,sk) = Vtest(:,i).*Vtest(:,j);
%                 sk = sk + 1;
%             end
%         end
%         b = sum(Vmass'*dx,2);
        
        % target space
        sk = 1;
        for i = 1:size(Vr,2)
            for j = i:size(Vtest,2)
                Vmass(:,sk) = Vr(:,i).*Vtest(:,j);
                sk = sk + 1;
            end
        end
        b = sum(Vmass'*dx,2);
        
%         % target space
%         sk = 1;
%         for i = 1:size(Vr,2)
%             for j = i:size(Vr,2)
%                 Vmass(:,sk) = Vr(:,i).*Vr(:,j);
%                 sk = sk + 1;
%             end
%         end
%         b = sum(Vmass'*dx,2);
        
        id = get_empirical_cubature(Vmass,b,.1*tol,5*Nmodes);
%         id = get_empirical_cubature(Vmass,b,tol,10*Nmodes);
        wr = Vmass(id,:)'\b;

        % make new mass, projection, interp matrices
        Mr = (Vr(id,:)'*diag(wr)*Vr(id,:));
        invM = inv(Mr);
        Pr = Mr\(Vr(id,:)'*diag(wr));
        Vr = Vr(id,:); 
        
        Mtest = Vtest(id,:)'*diag(wr)*Vtest(id,:);
        Ptest = Mtest\(Vtest(id,:)'*diag(wr));        
        S = Ptest'*(Vtest'*Sfull*Vtest)*Ptest;
        KS = Ptest'*(Vtest'*Kfull*Vtest)*Ptest;
    else
        S = Vtest*(Vtest'*Sfull*Vtest)*Vtest';
        invM = Vr'*Vr*(1/dx);
        Pr = invM*(Vr'*dx);        
    end

else
    
    invM = eye(K)/dx;
    Vr = eye(K);
    Pr = eye(K);
    Vrp = eye(K);
    
end


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

figure(1)
sk = 1;
for i = 1:Nsteps
    for INTRK = 1:5              
        
        rhoq = Vr*rho;
        mq = Vr*m;
        Eq = Vr*E; 
        
        v1 = Pr*V1(rhoq,mq,Eq);
        v2 = Pr*V2(rhoq,mq,Eq);
        v3 = Pr*V3(rhoq,mq,Eq);
        
        v1q = Vr*v1;
        v2q = Vr*v2;
        v3q = Vr*v3;
                      
        rhoq = U1(v1q,v2q,v3q);
        mq = U2(v1q,v2q,v3q);
        Eq = U3(v1q,v2q,v3q);
        uq = mq./rhoq;
                
        [rhox rhoy] = meshgrid(rhoq);
        [ux uy] = meshgrid(uq);
        [Ex Ey] = meshgrid(Eq);
        
        rhs1 = Vr'*sum(S.*fS1(rhox,rhoy,ux,uy,Ex,Ey),2);
        rhs2 = Vr'*sum(S.*fS2(rhox,rhoy,ux,uy,Ex,Ey),2);
        rhs3 = Vr'*sum(S.*fS3(rhox,rhoy,ux,uy,Ex,Ey),2);       
        
        if INTRK==5
           rhstest(i) = full(sum(sum(v1.*rhs1+v2.*rhs2+v3.*rhs3)));
        end                
        
        % add dissipation + invert mass
        tau = .5*dx;
        rhs1 = -invM*(rhs1 + tau*Vr'*KS*rhoq);
        rhs2 = -invM*(rhs2 + tau*Vr'*KS*mq);
        rhs3 = -invM*(rhs3 + tau*Vr'*KS*Eq);

        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        rho = rho  + rk4b(INTRK)*res1;
        m   = m  + rk4b(INTRK)*res2;
        E   = E  + rk4b(INTRK)*res3;
    end      
    
    err(i) = sqrt(sum(sum(dx*(Usnap(:,i+1)-[Vrp*rho;Vrp*m;Vrp*E]).^2)));
    
    if mod(i,10) == 0
        plot(x(id),Vr*rho,'o')
        hold on
        plot(x,Usnap(1:K,i+1),'-')
        hold off
        %         axis([-1,1,.5,3.5])
        
        title(sprintf('time = %f, step %d / %d, rhstest = %g, err = %g\n',i*dt,i,Nsteps,rhstest(i),err(i)));
        drawnow
    end
end

figure(2)
hold on
semilogy((1:Nsteps)*dt,err,'--','linewidth',2)

sqrt(sum(sum(dx*(Usnap(:,end)-[Vrp*rho;Vrp*m;Vrp*E]).^2)))
