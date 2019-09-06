% clear

K = 2500;
FinalTime = .75;
dx = 1/K;


Nmodes = 25;
tau = .5*dx;
snapshot_file = 'Usnap_euler_wall';
% snapshot_file = ['Usnap_euler_wall_K' num2str(K) '_taup5dx'];

xv = linspace(0,1,K+1)';
x = .5*(xv(1:end-1)+xv(2:end));

% tau = .1*dx^2;

e = ones(K-1,1);
S = diag(e,1)-diag(e,-1);
S(1,:) = 0; S(:,1) = 0;
S(end,:) = 0; S(:,end) = 0;
% S(1,end) = -1; S(end,1) = 1;
S(1,2) = 1; S(2,1) = -1;
S(K-1,K) = 1; S(K,K-1) = -1;

S = sparse(S);
KS = sparse(2*eye(K) - diag(e,1)-diag(e,-1));
KS(1,[1 2]) = [1 -1];
KS(K,[K-1 K]) = [-1 1];
KS = KS/dx;

B = spdiag([-1; zeros(K-2,1); 1]);

% for viscous hyper-rreduc
e = ones(K,1);
D = diag(e) - diag(e(2:end),1);
D = sparse(D(1:end-1,:));
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

% fS1 = @(rhoL,rhoR,uL,uR,EL,ER) avg(rhoL.*uL,rhoR.*uR);
% fS2 = @(rhoL,rhoR,uL,uR,EL,ER) avg(pfun(rhoL,uL,EL) + rhoL.*uL.^2,pfun(rhoR,uR,ER) + rhoR.*uR.^2);
% fS3 = @(rhoL,rhoR,uL,uR,EL,ER) avg(uL.*(EL+pfun(rhoL,uL,EL)),uR.*(ER+pfun(rhoR,uR,ER)));

%% make ROM

load(snapshot_file) 
    
U0 = reshape(Usnap(:,1),K,3);

Nsample = 1; %ceil(size(Usnap,2)/500); % downsample
Us1 = Usnap((1:K),1:Nsample:end);
Us2 = Usnap((1:K)+K,1:Nsample:end);
Us3 = Usnap((1:K)+2*K,1:Nsample:end);

% add entropy variables to snapshots
Us = [Us1 Us2 Us3 V1(Us1,Us2,Us3) V2(Us1,Us2,Us3) V3(Us1,Us2,Us3)];
tic;
[Vr,Sr,~] = rsvd(Us,5*Nmodes);
% [Vr,Sr,~] = svd(Us,0);
fprintf('Time for generating Vr = %f\n',toc)

Vr = Vr(:,1:Nmodes);
Vrfull = Vr;

sig = diag(Sr);
tol = sqrt(sum(sig(Nmodes+1:end).^2)/sum(sig.^2))

% add range + constants to snapshots
tic;[Vtest Stest, ~] = svd([ones(size(x)) Vr S*Vr],0);
fprintf('Time for generating Vtest = %f\n',toc)
sigt = diag(Stest);
Vtest = Vtest(:,sigt > 1e-12);

Sfull = S;
Kfull = KS;

rho = Vr'*U0(:,1);
m = Vr'*U0(:,2);
E = Vr'*U0(:,3);

rho0 = rho;
m0 = m;
E0 = E;

% % set tol based on init condition
% tol = sqrt(sum(sum(dx*(U0-Vr*[rho m E]).^2)))/sqrt(sum(sum(dx*U0.^2)))

%% empirical cubature
if 1        
    
    Vtest1 = Vrfull;
    Vtest2 = Vrfull; %
%     Vtest2 = Vtest;
    % target space
    tic;
    Vmass = zeros(size(Vr,1),size(Vr,2)*(size(Vr,2)+1)/2);
    sk = 1;
    for i = 1:size(Vtest1,2)
        for j = i:size(Vtest2,2)
            Vmass(:,sk) = Vtest1(:,i).*Vtest2(:,j);            
            sk = sk + 1;
        end
    end
    fprintf('Time for generating Vmass = %f\n',toc)
    
    % reduce target space
    tic;
    [Vmass,Smass,~] = rsvd(Vmass,25*Nmodes);
    %     [Vmass,Smass,~] = svd(Vmass,0);
    toc
    smass = diag(Smass);
    smass_energy = sqrt(1 - (cumsum(smass.^2)./sum(smass.^2)));
    Vmass = Vmass(:,smass_energy > tol);
    %     Vmass = Vmass(:,smass_energy > 1e-12);
    
    %     keyboard
    
    w0 = ones(K,1)*dx;
    [wr id] = get_empirical_cubature(Vmass,w0,tol,10*Nmodes);
    %     [wr id] = get_empirical_cubature(Vmass,ones(K,1)*dx,1e-12,10*Nmodes);
    
    if 0
        % gappy POD approach
        condV = cond(Vtest(id,:));
        if condV > 1e10
            fprintf('initial cond of Vtest(id,:) = %g\n',condV)
            sk = 0;
            while condV > 1e10
                ids = setdiff(1:K,id);
                for i = 1:length(ids)
                    idi = [id; i];
                    condV(i) = cond(Vtest(idi,:)'*Vtest(idi,:));
                    if mod(i,500)==0
                        fprintf('on i = %d out of %d\n',i,length(ids))
                    end
                end
                [condV,i] = min(condV);
                id = [id;ids(i)];
                sk = sk + 1;               
                fprintf('adding %d points, condV = %g\n',sk,condV);
            end            
            keyboard
        end
        Ptest = pinv(Vtest(id,:)); 
    
    else
        
        % L2 projection approach
        Mtest = Vtest(id,:)'*diag(wr)*Vtest(id,:);
        [Vx Dx] = eig(Mtest);
        d = diag(Dx);
        [d,p] = sort(d,'ascend');
        Z = Vtest*Vx(:,p(d < 1e-12));        
                
        if cond(Mtest) > 1e10
            fprintf('initial cond of Mtest %g\n',cond(Mtest))
            Zmass = zeros(size(Z,1),size(Z,2)*(size(Z,2)+1)/2);
            sk = 1;
            for i = 1:size(Z,2)
                for j = 1:size(Z,2)
                    Zmass(:,sk) = Z(:,i).*Z(:,j);
                    sk = sk + 1;
                end
            end
            [Zmass,SZ,~] = svd(Zmass,0);
            Zmass = Zmass(:,diag(SZ) > 1e-12);
%             Zmass = Zmass(:,diag(SZ) > tol);
            fprintf('\n getting additional points for stability\n')
            [wz idz] = get_empirical_cubature(Zmass,w0,tol,size(Zmass,2));
            
            zscale = 1e-2;
            idt = unique([id; idz]);
            J = [Vmass zscale*Zmass];
            b = sum(J',2);
            w_unscaled = lsqlin(J(idt,:)',b,[],[],[] ,[] ,zeros(size(idt)),inf(size(idt)),J(idt,:)'\b); % upper bound = total vol
            w = w_unscaled*sum(w0)/K;
            
            Mtest = Vtest(idt,:)'*diag(w)*Vtest(idt,:);
%             [Vt Dt] = eig(Mtest);
%             lam = diag(Dt);
%             [minlam,idx] = min(lam);
            
            if cond(Mtest) > 1e10
                keyboard
            end
            
            id = idt;
            wr = w;
            
            fprintf('Added %d new points, new cond of Mtest = %g\n',length(idz),cond(Mtest))
            Mtest = Vtest(id,:)'*diag(wr)*Vtest(id,:);
            Ptest = Mtest\(Vtest(id,:)'*diag(wr));
        end
        
        
    end
    
    % make new mass, projection, interp matrices
    Mr = (Vr(id,:)'*diag(wr)*Vr(id,:));
    invM = inv(Mr);
    Pr = Mr\(Vr(id,:)'*diag(wr));
        
    Vr = Vr(id,:);
    
    KS = Vr'*Pr'*(Vrfull'*Kfull*Vrfull)*Pr;
        
    S = Ptest'*(Vtest'*Sfull*Vtest)*Ptest;        
    
    % provably stable hyper-reduction for viscous terms
    tic;[VD, SD, ~] = rsvd(D*Vrfull,25*Nmodes);
    smass = diag(SD);
    smass_energy = sqrt(1 - (cumsum(smass.^2)./sum(smass.^2)));
    VD = VD(:,smass_energy > tol);
    fprintf('Computed rsvd of D*V in %g sec\n',toc);
    
    fprintf('Forming VmassD\n')
    tic
    VmassD = zeros(size(VD,1),size(VD,2)*(size(VD,2)+1)/2);
    sk = 1;
    for i = 1:size(VD,2)
        for j = 1:size(VD,2)
            VmassD(:,sk) = VD(:,i).*VD(:,j);
            sk = sk + 1;
        end
    end    
    [Vmass,Smass,~] = rsvd(VmassD,25*Nmodes);
    toc
    smass = diag(Smass);
    smass_energy = sqrt(1 - (cumsum(smass.^2)./sum(smass.^2)));
    Vmass = Vmass(:,smass_energy > tol);
    [wDr idD] = get_empirical_cubature(Vmass,ones(size(D,1),1),tol,15*Nmodes); % cheaper
    
    %         wDr = ones(size(D,1),1);
    %         idD = 1:size(D,1);
    
    DV = D*Vrfull;
    DV = DV(idD,:)/sqrt(dx);  % DV'*DV is about the same as Vrp'*Kfull*Vrp
    VrD = .5*(Vrfull(idD,:)+Vrfull(idD+1,:)); % eval @ midpoints
    
    
    [i j] = find(D(idD,:));
    %         norm(D(idD,unique(j))*Vrp(unique(j),:) - D(idD,:)*Vrp,'fro')
    Dsparse = D(idD,unique(j))/sqrt(dx);
    VrDsparse = Vrfull(unique(j),:);
    
    % lax friedrichs
    for i = 1:size(Dsparse,1)
        ids = find(Dsparse(i,:));
        idDL(i) = ids(1);
        idDR(i) = ids(2);
    end
%     keyboard
    
else
    
    invM = eye(K)/dx;
    Pr = eye(K);
    Vr = eye(K);
    KS = Kfull;
    S = Sfull;
end

Vf = Vrfull([1 K],:);
Vftest = Vtest([1 K],:);
Eftest = Vftest*Ptest;
B = diag([-1;1]);

QN = [S Eftest'*B;
    -B*Eftest zeros(2)];

VN = [Vr; Vf];

%% 

rho = rho0;
m = m0;
E = E0;

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
res1 = zeros(size(rho));
res2 = zeros(size(rho));
res3 = zeros(size(rho));

figure
for i = 1:Nsteps
    for INTRK = 1:5                        
        
        rhoq = Vr*rho;
        mq = Vr*m;
        Eq = Vr*E;
        v1N = Pr*V1(rhoq,mq,Eq);
        v2N = Pr*V2(rhoq,mq,Eq);
        v3N = Pr*V3(rhoq,mq,Eq);
        v1 = VN*v1N;
        v2 = VN*v2N;
        v3 = VN*v3N;
        
        rhoN = U1(v1,v2,v3);
        mN = U2(v1,v2,v3);
        EN = U3(v1,v2,v3);
        uN = mN./rhoN;
        betaN = beta(rhoN,uN,EN);
        
        [rhox rhoy] = meshgrid(rhoN);
        [ux uy] = meshgrid(uN);
        [bx by] = meshgrid(betaN);                
        
        rhoq = rhoN(1:end-2);
        mq = mN(1:end-2);
        Eq = EN(1:end-2);
        
        rhof = rhoN(end-1:end);
        uf = uN(end-1:end);
        Ef = EN(end-1:end);
%         uf = mf./rhof;        
        
        % fluxes + BCs
        rhoL = rhof;
        rhoR = rhof;
        EL = Ef;
        ER = Ef;
        uL = [-uf(1);uf(2)];
        uR = [uf(1);-uf(2)];                
        
        % pavg = @(rhoL,rhoR,uL,uR,EL,ER) avg(rhoL,rhoR)./(2*avg(beta(rhoL,uL,EL),beta(rhoR,uR,ER)));
        % fS1 = @(rhoL,rhoR,uL,uR,EL,ER) logmean(rhoL,rhoR).*avg(uL,uR);
        % fS2 = @(rhoL,rhoR,uL,uR,EL,ER) pavg(rhoL,rhoR,uL,uR,EL,ER) + avg(uL,uR).*fS1(rhoL,rhoR,uL,uR,EL,ER);
        % fS3 = @(rhoL,rhoR,uL,uR,EL,ER) fS1(rhoL,rhoR,uL,uR,EL,ER)...
        %     .*(1./(2*(gamma-1).*logmean(beta(rhoL,uL,EL),beta(rhoR,uR,ER))) - .5*avg(uL.^2,uR.^2)) ...
        %     + avg(uL,uR).*fS2(rhoL,rhoR,uL,uR,EL,ER);
        
        rholog = logmean(rhox,rhoy);
        betalog = logmean(bx,by);
        uavg = avg(ux,uy);
        pa = avg(rhox,rhoy)./(2*avg(bx,by));
        uavg2 = 2*uavg.^2-avg(ux.^2,uy.^2);
        Ea = rholog.*(1./(2*(gamma-1)*betalog)+uavg2/2);
        F1 = rholog.*uavg;
        F2 = pa + uavg.*F1;
        F3 = (Ea+pa).*uavg;
        
%         F1 = fS1(rhox,rhoy,ux,uy,Ex,Ey);
%         F2 = fS2(rhox,rhoy,ux,uy,Ex,Ey);
%         F3 = fS3(rhox,rhoy,ux,uy,Ex,Ey);
        rhs1 = VN'*sum(QN.*F1,2) + Vf'*([-1;1].*fS1(rhoL,rhoR,uL,uR,EL,ER));
        rhs2 = VN'*sum(QN.*F2,2) + Vf'*([-1;1].*fS2(rhoL,rhoR,uL,uR,EL,ER));
        rhs3 = VN'*sum(QN.*F3,2) + Vf'*([-1;1].*fS3(rhoL,rhoR,uL,uR,EL,ER));  
        
        if (INTRK==5)
            rhstest(i) = dx*sum(sum(v1N'*rhs1 + v2N'*rhs2 + v3N'*rhs3));
        end
        
        % Lax-Friedrichs penalty
        tau_lf = 1/2;
        unorm2 = (uf.^2);
        pf = (gamma-1)*(Ef - .5*rhof.*unorm2);
        cvel = sqrt(gamma*pf./rhof);
        lam = sqrt(unorm2)+cvel;
        LFc = abs(lam);        
        % rhoP = rhoM, EP = EM so LF only activates 2nd penalty = mP= -mM
        LF2 = tau_lf*LFc.*(rhof.*(-2*uf));         
        rhs2 = rhs2 - Vf'*LF2;
        
        opt = 2;
        if opt==3
            d1 = KS*rhoq;
            d2 = KS*mq;
            d3 = KS*Eq;
        elseif opt==2
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
        elseif opt==1
            v1 = VrDsparse*v1N;
            v2 = VrDsparse*v2N;
            v3 = VrDsparse*v3N;
            
            rhoq = U1(v1,v2,v3);
            mq   = U2(v1,v2,v3);
            Eq   = U3(v1,v2,v3);
            d1 = VrDsparse'*Dsparse'*(wDr.*(Dsparse*rhoq));
            d2 = VrDsparse'*Dsparse'*(wDr.*(Dsparse*mq));
            d3 = VrDsparse'*Dsparse'*(wDr.*(Dsparse*Eq));
        end
        
        % add dissipation + invert mass        
        rhs1 = -invM*(rhs1 + tau*d1);
        rhs2 = -invM*(rhs2 + tau*d2);
        rhs3 = -invM*(rhs3 + tau*d3);
        
        dtest(i) = v1N'*d1 + v2N'*d2 + v3N'*d3;

        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        rho = rho  + rk4b(INTRK)*res1;
        m   = m  + rk4b(INTRK)*res2;
        E   = E  + rk4b(INTRK)*res3;
    end
    
    err(i) = sqrt(sum(sum(dx*(reshape(Usnap(:,i),K,3)-Vrfull*[rho m E]).^2)));
    
    Sq = -(Vrfull*rho).*s(Vrfull*rho,Vrfull*m,Vrfull*E);
    rom_entropy(i) = sum(sum(dx.*Sq));
    
    if mod(i,10) == 0 || i==Nsteps
        plot(x(id),Vr*rho,'ro','linewidth',2)
        hold on
        plot(x,Usnap(1:K,i),'b-','linewidth',2)
        hold off
        title(sprintf('time = %f, step %d / %d, err = %g, rhstest = %g, dissip = %g',dt*i,i,Nsteps,err(i),rhstest(i),dtest(i)));
        drawnow
    end
end

axis([0 1 1.95 2.4]); grid on; title(''); set(gca,'fontsize',16);
% print(gcf,'-dpng',sprintf('~/Desktop/bbwadg/docs/ESDG_ROM/figs/euler1Dwall_t%d_%dmodes.png',round(100*FinalTime),Nmodes))

% plot(x,Vrfull*rho-Usnap(1:K,i))
sqrt(sum(sum(dx*(reshape(Usnap(:,Nsteps),K,3)-Vrfull*[rho m E]).^2)))/sqrt(sum(sum(dx*(reshape(Usnap(:,Nsteps),K,3)).^2)))
% rho(end-10:end) = 0;




figure(3)
semilogy(dt*(1:Nsteps),dtest,'o--','markersize',10)
hold on
grid on

fontsz = 16;
xlabel('Time','fontsize',fontsz)
ylabel({'Discrete entropy dissipation'},'fontsize',fontsz);
% legend show
% legend location southwest
set(gca,'fontsize',fontsz)


return
%%

figure(3)
plot(x(id),Vr*rho,'ro','linewidth',2)
hold on
plot(x,Usnap(1:K,i),'b-','linewidth',2)
set(gca,'xticklabel','','yticklabel','')
grid on

% print(gcf,'-dpng','~/Desktop/proposals/CAREER/figs/euler1Dwall_t75_125modes.png')
%%
figure(4)
semilogy(dt*(1:Nsteps),fom_entropy+0*fom_entropy(end),'b-','linewidth',2,'DisplayName','Full order model')
% semilogy(dt*(1:Nsteps),abs(fom_entropy-rom_entropy),'b-','linewidth',2,'DisplayName','Full order model')
hold on
semilogy(dt*(1:Nsteps),rom_entropy+0*fom_entropy(end),'r--','linewidth',2,'DisplayName',sprintf('Reduced model (%d modes)',Nmodes))
grid on

fontsz = 16;
xlabel('Time','fontsize',fontsz)
ylabel({'Average entropy'},'fontsize',fontsz);
legend show
% legend location southwest
set(gca,'fontsize',fontsz)

% print(gcf,'-dpng','~/Desktop/proposals/CAREER/figs/euler1Dwall_entropy_75modes.png')
