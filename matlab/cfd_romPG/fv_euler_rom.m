clear
K = 1000;
FinalTime = .7;
xv = linspace(0,1,K+1)';
x = .5*(xv(1:end-1)+xv(2:end));
dx = 1/K;

tau = .5*dx;
% tau = 1e-3;
% tau = .25*dx;

e = ones(K-1,1);
S = diag(e,1)-diag(e,-1);
S(1,:) = 0; S(:,1) = 0;
S(end,:) = 0; S(:,end) = 0;
S(1,end) = -1; S(end,1) = 1;
S(1,2) = 1; S(2,1) = -1;
S(K-1,K) = 1; S(K,K-1) = -1;

S = sparse(S);
% KS = sparse(2*eye(K) - abs(S))/dx^2;

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

% entropy potentials
psi = @(rho,rhou,E) (gamma-1)*rhou;


%% construct ROM

if 1
    load Usnap_euler
%     load Fsnap_euler
    
    U0 = reshape(Usnap(:,1),K,3);
    %  rho = Vr'*rho;
    %  m = Vr'*m;
    %  E = Vr'*E;    

    Nsample = 2; %ceil(size(Usnap,2)/500); % downsample
    Us1 = Usnap((1:K),1:Nsample:end);
    Us2 = Usnap((1:K)+K,1:Nsample:end);
    Us3 = Usnap((1:K)+2*K,1:Nsample:end);
    
    % add entropy variables to snapshots
    Us = [Us1 Us2 Us3 V1(Us1,Us2,Us3) V2(Us1,Us2,Us3) V3(Us1,Us2,Us3)];
%     Us = [Us1 Us2 Us3];
    [Vr,Sr,~] = svd(Us,0);
    
    Nmodes = 25;
    Vr = Vr(:,1:Nmodes);
    Vrp = Vr;
    
    sig = diag(Sr);
    tol = sqrt(sum(sig(Nmodes+1:end).^2)/sum(sig.^2))
    
    % add range + constants to snapshots
    [Vtest Stest, ~] = svd([Vr S*Vr]);
    sigt = diag(Stest);        
    Vtest = Vtest(:,sigt > 1e-10);
    Vtest = orth([ones(size(x)) Vtest]); % ensure 1
    [Vrange, Srange, ~] = svd(S*Vr,0);
    Vrange = Vrange(:,diag(Srange)>1e-10);
        
    Sfull = S;
    Kfull = KS;
       
    rho = Vr'*U0(:,1);
    m = Vr'*U0(:,2);
    E = Vr'*U0(:,3);   
    
%     % set tol based on init condition
%     tol = sqrt(sum(sum(dx*(U0-Vrp*[rho m E]).^2)))/sqrt(sum(sum(dx*U0.^2)))
    
    % empirical cubature
    if 1        
        % target space
        Vtest1 = Vr;         
        Vtest2 = Vr; 
%         Vtest2 = Vtest; 
        sk = 1;        
        for i = 1:size(Vtest1,2)
            for j = 1:size(Vtest2,2)
                Vmass(:,sk) = Vtest1(:,i).*Vtest2(:,j);
                sk = sk + 1;
            end
        end        
        
        % reduce target space
        [Vmass,Smass,~] = svd(Vmass,0);
        smass = diag(Smass);
        smass_energy = sqrt(1 - (cumsum(smass.^2)./sum(smass.^2)));                
        if nnz(smass_energy > tol) > size(Vtest,2)
            Vmass = Vmass(:,smass_energy > tol);
        else
            Vmass = Vmass(:,1:size(Vtest,2));
        end
%         Vmass = Vmass(:,smass_energy > 0);
                        

        [wr id] = get_empirical_cubature(Vmass,ones(size(x(:)))*dx,tol,15*Nmodes); % cheaper
%  [wr id] = get_EQP_nodes(Vmass,ones(size(x(:)))*dx,tol);
        
        % make new mass, projection, interp matrices
        Mr = (Vr(id,:)'*diag(wr)*Vr(id,:));
        invM = inv(Mr);
        Pr = Mr\(Vr(id,:)'*diag(wr));
        Vr = Vr(id,:); 
        
%         Ptest = (Vtest(id,:)'*diag(wr)*Vtest(id,:))\(Vtest(id,:)'*diag(wr));        
%         Ptest = (Vtest(id,:)'*Vtest(id,:))\(Vtest(id,:)');     
        Vr1 = [Vrp ones(size(Vrp,1),1)]; % augment with 1
        R = Vtest\Vr1;
        Pr1 = (Vr1(id,:)'*diag(wr)*Vr1(id,:)) \ (Vr1(id,:)'*diag(wr));
        Ptest = R*Pr1;
        S = Ptest'*(Vtest'*Sfull*Vtest)*Ptest;
        KS = Ptest'*(Vtest'*Kfull*Vtest)*Ptest;
        VrKS = Vrp'*Kfull*Vrp*Pr;
        
        % provably stable hyper-reduction for viscous terms        
        [VD, SD, ~] = svd(D*Vrp,0);
        VD = VD(:,diag(SD) > tol);%1e-12);
        sk = 1;        
        for i = 1:size(VD,2)
            for j = 1:size(VD,2)
                VmassD(:,sk) = VD(:,i).*VD(:,j);
                sk = sk + 1;
            end
        end        
        [Vmass,SmassD,~] = svd(VmassD,0);
        smass = diag(SmassD);
        smass_energy = sqrt(1 - (cumsum(smass.^2)./sum(smass.^2)));                
        Vmass = Vmass(:,smass_energy > tol);
        [wDr idD] = get_empirical_cubature(Vmass,ones(size(D,1),1),tol,15*Nmodes); % cheaper
        
%         wDr = ones(size(D,1),1);
%         idD = 1:size(D,1);
        
        DV = D*Vrp;
        DV = DV(idD,:)/sqrt(dx);  % DV'*DV is about the same as Vrp'*Kfull*Vrp
        VrD = .5*(Vrp(idD,:)+Vrp(idD+1,:)); % eval @ midpoints
        
        [i j] = find(D(idD,:));
%         norm(D(idD,unique(j))*Vrp(unique(j),:) - D(idD,:)*Vrp,'fro')        
        Dsparse = D(idD,unique(j))/sqrt(dx);
        VrDsparse = Vrp(unique(j),:);
        
        % lax friedrichs
        for i = 1:size(Dsparse,1)
            ids = find(Dsparse(i,:));
            idDL(i) = ids(1);
            idDR(i) = ids(2);
        end
            
%         VrD = Vrp(idD,:);
        
%         norm(Vrp'*Kfull*Vrp - DV'*DV,'fro')
        
%         % exact mass matrix
%         invM = eye(size(Vr,2))/dx;
%         Pr = invM*Vr'*diag(wr);
% keyboard
        
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
        
        v1N = Pr*V1(rhoq,mq,Eq);
        v2N = Pr*V2(rhoq,mq,Eq);
        v3N = Pr*V3(rhoq,mq,Eq);
        
        v1q = Vr*v1N;
        v2q = Vr*v2N;
        v3q = Vr*v3N;
                      
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
%             keyboard
           rhstest(i) = full(dx*sum(sum(v1N.*rhs1+v2N.*rhs2+v3N.*rhs3)));
        end                
        
        % add dissipation + invert mass
        
        opt=2;
        if opt==0
            d1 = VrKS*rhoq;
            d2 = VrKS*mq;
            d3 = VrKS*Eq;
        elseif opt==1        
            % D*V*vN
            v1D = DV*v1N;
            v2D = DV*v2N;
            v3D = DV*v3N;
            
            % evaluate entropy var at hyper-reduced points            
            v1 = VrD*v1N;
            v2 = VrD*v2N;
            v3 = VrD*v3N;
            scale = rhoeV(v1,v2,v3)./((gamma-1)*v3);
            v23 = v2.^2./(2*v3);
            d11 = scale.*(-v3.^2);
            d12 = scale.*(v2.*v3);
            d13 = scale.*(v3.*(1-v23));
            d22 = scale.*((gamma-1)*v3-v2.^2);
            d23 = scale.*(v2.*(v23-gamma));
            d33 = scale.*(-((v2.^2./(2*v3)).^2-2*gamma*(v2.^2./(2*v3)) + gamma));
            
            d1 = DV'*(wDr.*(d11.*v1D + d12.*v2D + d13.*v3D));
            d2 = DV'*(wDr.*(d12.*v1D + d22.*v2D + d23.*v3D));
            d3 = DV'*(wDr.*(d13.*v1D + d23.*v2D + d33.*v3D)); 
            
        elseif opt==2
            
            v1 = VrDsparse*v1N;
            v2 = VrDsparse*v2N;
            v3 = VrDsparse*v3N;
            
            rhoq = U1(v1,v2,v3);
            mq   = U2(v1,v2,v3);
            Eq   = U3(v1,v2,v3);
            d1 = VrDsparse'*Dsparse'*(wDr.*(Dsparse*rhoq));
            d2 = VrDsparse'*Dsparse'*(wDr.*(Dsparse*mq));
            d3 = VrDsparse'*Dsparse'*(wDr.*(Dsparse*Eq));
        
        elseif opt==3 % lax friedrichs penalization
            
            v1 = VrDsparse*v1N;
            v2 = VrDsparse*v2N;
            v3 = VrDsparse*v3N;
            
            rhoq = U1(v1,v2,v3);
            mq   = U2(v1,v2,v3);
            Eq   = U3(v1,v2,v3);
            
            rhoL = rhoq(idDL);
            rhoR = rhoq(idDL+1);
            uL = mq(idDL)./rhoL;
            uR = mq(idDL+1)./rhoR;
            EL = Eq(idDL);
            ER = Eq(idDL+1);
            
            % tau is scaled by dx since LF flux isn't scaled
            Lfc = llf(rhoL,rhoR,uL,uR,EL,ER);
            d1 = VrDsparse'*Dsparse'*(wDr.*Lfc.*(Dsparse*rhoq));
            d2 = VrDsparse'*Dsparse'*(wDr.*Lfc.*(Dsparse*mq));
            d3 = VrDsparse'*Dsparse'*(wDr.*Lfc.*(Dsparse*Eq));
        end
        
        if (INTRK==5)
            dtest(i) = v1N'*d1 + v2N'*d2 + v3N'*d3;
        end
        
        rhs1 = -invM*(rhs1 + tau*d1);
        rhs2 = -invM*(rhs2 + tau*d2);
        rhs3 = -invM*(rhs3 + tau*d3);

        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        rho = rho  + rk4b(INTRK)*res1;
        m   = m  + rk4b(INTRK)*res2;
        E   = E  + rk4b(INTRK)*res3;
    end      
    
    err(i) = sqrt(sum(sum(dx*(Usnap(:,i+1)-[Vrp*rho;Vrp*m;Vrp*E]).^2)));
    Sq = -rho.*s(rho,m,E);
    entropy(i) = sum(sum(dx.*Sq));
    
    if mod(i,5) == 0
        plot(x(id),Vr*rho,'o')
        hold on
        plot(x,Usnap(1:K,i+1),'-')
        hold off
        %         axis([-1,1,.5,3.5])
        
        title(sprintf('time = %f, step %d / %d, rhstest = %g, dtest = %g, err = %g\n',i*dt,i,Nsteps,rhstest(i),dtest(i),err(i)));
        drawnow
    end
end

figure(2)
semilogy((1:Nsteps)*dt,err,'--','linewidth',2)
hold on

sqrt(sum(sum(dx*(Usnap(:,Nsteps+1)-[Vrp*rho;Vrp*m;Vrp*E]).^2)))

figure(3)
semilogy(dt*(1:Nsteps),rhstest,'o--')
hold on

% figure(3)
% hold on
% semilogy((1:Nsteps)*dt,abs(entropy-entropy(1)),'--','linewidth',2)
