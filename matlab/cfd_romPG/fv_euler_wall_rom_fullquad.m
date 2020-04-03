% clear
K = 2500;
FinalTime = .75;

Nmodes = 200;
tau = .5*dx;
snapshot_file = 'Usnap_euler_wall';
% snapshot_file = ['Usnap_euler_wall_K' num2str(K) '_taup5dx'];
% snapshot_file = ['Usnap_euler_wall_K' num2str(K) '_tau0'];


xv = linspace(0,1,K+1)';
x = .5*(xv(1:end-1)+xv(2:end));
dx = 1/K;

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

%% set initial cond

opt = 0;

if opt==0    
    % sine solution
    t = 0;    
    a = 1;
    rhoex = @(x) 2 + a*exp(-100*(x-.5).^2);    
    uex = @(x) zeros(size(x));
    pex = @(x) rhoex(x).^gamma;    
    
elseif opt==1
    % afkam, hesthaven
    
    rhoex = @(x) .5 + .2*cos(2*pi*x);
    uex = @(x) 1.5 + 0*x;
    pex = @(x) .5 + .2*sin(2*pi*x);
end

mex = @(x) rhoex(x).*uex(x);
Eex = @(x) pex(x)./(gamma-1) + .5*rhoex(x).*uex(x).^2;

rho = rhoex(x);
m = mex(x);
E = Eex(x);

% plot(x,E,'o')
% return

%% make ROM

load(snapshot_file);
    
U0 = reshape(Usnap(:,1),K,3);

Nsample = 2; %ceil(size(Usnap,2)/500); % downsample
Us1 = Usnap((1:K),1:Nsample:end);
Us2 = Usnap((1:K)+K,1:Nsample:end);
Us3 = Usnap((1:K)+2*K,1:Nsample:end);

% add entropy variables to snapshots
Us = [Us1 Us2 Us3 V1(Us1,Us2,Us3) V2(Us1,Us2,Us3) V3(Us1,Us2,Us3)];
tic;[Vr,Sr,~] = svd(Us,0);
fprintf('Time for generating Vr = %f\n',toc)

Vr = Vr(:,1:Nmodes);
Vrfull = Vr;

sig = diag(Sr);
tol = sqrt(sum(sig(Nmodes+1:end).^2)/sum(sig.^2))

% add range + constants to snapshots
tic;[Vtest Stest, ~] = svd([ones(size(x)) Vr S*Vr],0);
fprintf('Time for generating Vtest = %f\n',toc)
sigt = diag(Stest);
Vtest = Vtest(:,sigt > 1e-13);

% [Vrange,Srange,~] = svd(S*Vr,0);
% Vrange = Vrange(:,diag(Srange) > 1e-13);

Sfull = S;
Kfull = KS;

rho = Vr'*U0(:,1);
m = Vr'*U0(:,2);
E = Vr'*U0(:,3);

% % set tol based on init condition
% tol = sqrt(sum(sum(dx*(U0-Vr*[rho m E]).^2)))/sqrt(sum(sum(dx*U0.^2)))

%% empirical cubature

if 1
    
%     Vtest1 = Vr;
%     Vtest2 = Vr; %
% %     Vtest2 = Vtest;
%     % target space
%     tic;Vmass = zeros(size(Vr,1),size(Vr,2)*(size(Vr,2)+1)/2);
%     sk = 1;
%     for i = 1:size(Vtest1,2)
%         for j = i:size(Vtest2,2)
%             Vmass(:,sk) = Vtest1(:,i).*Vtest2(:,j);            
%             sk = sk + 1;
%         end
%     end
%     fprintf('Time for generating Vmass = %f\n',toc)
%     
%     % reduce target space
%     tic;[Vmass,Smass,~] = rsvd(Vmass,50*Nmodes);toc
%     smass = diag(Smass);
%     smass_energy = sqrt(1 - (cumsum(smass.^2)./sum(smass.^2)));
% %     Vmass = Vmass(:,smass_energy > tol);
%     Vmass = Vmass(:,smass_energy > 1e-12);    
%     
%     [wr id] = get_empirical_cubature(Vmass,ones(K,1)*dx,tol,10*Nmodes);
% %     [wr id] = get_empirical_cubature(Vmass,ones(K,1)*dx,1e-12,10*Nmodes);

% full quad
wr = ones(K,1)*dx;
id = 1:K;
    
    % make new mass, projection, interp matrices
    Mr = (Vr(id,:)'*diag(wr)*Vr(id,:));
    invM = inv(Mr);
    Pr = Mr\(Vr(id,:)'*diag(wr));
    Vr = Vr(id,:);
    
    KS = Vr'*Pr'*(Vrfull'*Kfull*Vrfull)*Pr;
%     Mtest = Vtest(id,:)'*diag(wr)*Vtest(id,:);
%     Ptest = Mtest\(Vtest(id,:)'*diag(wr));
%     S = Ptest'*(Vtest'*Sfull*Vtest)*Ptest;        

else
    
    invM = eye(K)/dx;
    Pr = eye(K);
    Vr = eye(K);
    KS = Kfull;
    S = Sfull;
end

Vf = Vrfull([1 K],:);
% Vftest = Vtest([1 K],:);
% Eftest = Vftest*Ptest;
% B = diag([-1;1]);
% 
% QN = [S Eftest'*B;
%     -B*Eftest zeros(2)];
% VN = [Vr; Vf];


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
        v1 = Vr*v1N;
        v2 = Vr*v2N;
        v3 = Vr*v3N;
        
        rhoq = U1(v1,v2,v3);
        mq   = U2(v1,v2,v3);
        Eq   = U3(v1,v2,v3);
        uq = mq./rhoq;
        
        rhoL = [rhoq(1) rhoq(:)'];
        rhoR = [rhoq(:)' rhoq(end)];
        uL = [-uq(1) uq(:)'];
        uR = [uq(:)' -uq(end)];        
        EL = [Eq(1) Eq(:)'];
        ER = [Eq(:)' Eq(end)];                 
        
        r1 = Vr'*diff(fS1(rhoL,rhoR,uL,uR,EL,ER))';
        r2 = Vr'*diff(fS2(rhoL,rhoR,uL,uR,EL,ER))';
        r3 = Vr'*diff(fS3(rhoL,rhoR,uL,uR,EL,ER))';
        
        if INTRK==5
           rhstest(i) = full(sum(sum(v1N.*r1 + v2N.*r2 + v3N.*r3)));
        end 
        
        % Lax-Friedrichs penalty
        rhof = [rhoq(1);rhoq(end)];
        uf = [uq(1);uq(end)];
        Ef = [Eq(1);Eq(end)];
        tau_lf = 1/2;
        unorm2 = (uf.^2);
        pf = (gamma-1)*(Ef - .5*rhof.*unorm2);
        cvel = sqrt(gamma*pf./rhof);
        lam = sqrt(unorm2)+cvel;
        LFc = abs(lam);        
        % rhoP = rhoM, EP = EM so LF only activates 2nd penalty = mP= -mM
        LF1 = 0;
        LF2 = tau_lf*LFc.*(rhof.*(-2*uf)); 
        LF3 = 0;
        r2 = r2 - Vf'*LF2;
                                
        % add dissipation + invert mass
        rhs1 = -invM*(r1 + tau*KS*rhoq);
        rhs2 = -invM*(r2 + tau*KS*mq);
        rhs3 = -invM*(r3 + tau*KS*Eq);
        
        dissip = 0; %v1'*KS*rhoq+v2'*KS*mq + v3'*KS*Eq;

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
        title(sprintf('time = %f, step %d / %d, err = %g, rhstest = %g, dissip = %g',dt*i,i,Nsteps,err(i),rhstest(i),dissip));
        drawnow
    end
end

% plot(x,Vrfull*rho-Usnap(1:K,i))
sqrt(sum(sum(dx*(reshape(Usnap(:,Nsteps),K,3)-Vrfull*[rho m E]).^2)))/sqrt(sum(sum(dx*(reshape(Usnap(:,Nsteps),K,3)).^2)))
% rho(end-10:end) = 0;

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
