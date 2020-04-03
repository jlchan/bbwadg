clear
K = 300;
FinalTime = .1;
% FinalTime = .25;

xv = linspace(-1,1,K+1)';
x1D = .5*(xv(1:end-1)+xv(2:end));
x = x1D;
dx = 1/K;

CFL = .5;
Nmodes = 50;
tau = .1*dx;
tau = 2e-3;

% snapshot_file = 'Usnap_euler2d_wall';
snapshot_file = 'Usnap_euler2d_shock';

e = ones(K-1,1);
S = diag(e,1)-diag(e,-1);
S(1,:) = 0; S(:,1) = 0;
S(end,:) = 0; S(:,end) = 0;
% S(1,end) = -1; S(end,1) = 1;
S(1,2) = 1; S(2,1) = -1;
S(K-1,K) = 1; S(K,K-1) = -1;

S = sparse(S);
% KSx = sparse(2*eye(K) - abs(S))/dx;
% KS = sparse(2*eye(K) - diag(e,1) - diag(e,-1))/dx;
% KS(1,1) = KS(1,1)/2;
% KS(K,K) = KS(K,K)/2;

e = ones(K,1);
D = diag(e) - diag(e(2:end),1);
D = D(1:end-1,:);
KS = sparse(D'*D)/dx;


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

%% build 2D ops

% add range + constants to snapshots
Qxfull = kron(S,dx*speye(K));
Qyfull = kron(dx*speye(K),S);
% Kfull = kron(KS,dx*speye(K)) + kron(dx*speye(K),KS); % d2u/dx2 + d2u/dy2;

Dfullx = kron(D,speye(K));
Dfully = kron(speye(K),D);
% norm(Kfull-(Dfullx'*Dfullx + Dfully'*Dfully),'fro')
Kfull = (Dfullx'*Dfullx + Dfully'*Dfully);

%% make ROM

load(snapshot_file) 
% load Usnap_euler2d_wall_p5dxvisc
% load Usnap_euler2d_wall_dxvisc
% load Usnap_euler2d_wall_novisc

U0 = reshape(Usnap(:,1),K^2,4);
u1 = Usnap(1:K^2,:);
u2 = Usnap((1:K^2) + K^2,:);
u3 = Usnap((1:K^2) + 2*K^2,:);
u4 = Usnap((1:K^2) + 3*K^2,:);

[Vr, Sr, ~] = svd([u1 u2 u3 u4 ...
    V1(u1,u2,u3,u4) V2(u1,u2,u3,u4) V3(u1,u2,u3,u4) V4(u1,u2,u3,u4)],0); 
sig = diag(Sr); 

sigU = svd([u1 u2 u3 u4]);
semilogy(sigU,'bo','linewidth',2,'markersize',10)
hold on
semilogy(sig,'rx','linewidth',2,'markersize',10)
ylim([1e-12,1e5])
grid on
set(gca,'fontsize',15)
legend('Without enrichment','With enrichment')

keyboard
% figure(1)

% hold on
% grid on
% set(gca,'fontsize',16)


Vr = Vr(:,1:Nmodes);
Vrp = Vr;

tol = sum(sig(Nmodes+1:end).^2)/sum(sig.^2)

%% make test ops

DVrx = [Qxfull*Vr];
DVry = [Qyfull*Vr];

% reduce DVr size
tic; 
sDVr = svd(DVrx,0); toc
sDVr_energy = sqrt(1 - (cumsum(sDVr.^2)./sum(sDVr.^2)));
DVrx = DVrx(:,sDVr_energy > 1e-12); % this is important!  Don't do tol here.

sDVr = svd(DVry,0); toc
sDVr_energy = sqrt(1 - (cumsum(sDVr.^2)./sum(sDVr.^2)));
DVry = DVry(:,sDVr_energy > 1e-12); % this is important!  Don't do tol here.
disp('Done with DVr svd')

% construct Vtest by augmenting basis 
tic;
[Vtestx, Stest, ~] = svd([Vr DVrx ones(size(x(:)))],0);toc
sigt = diag(Stest);
Vtestx = Vtestx(:,sigt > 1e-12); 
% Vtestx = orth([ones(size(x(:))) Vtestx]); % ensure 1 in test space

[Vtesty, Stest, ~] = svd([Vr DVry ones(size(x(:)))],0);toc
sigt = diag(Stest);
Vtesty = Vtesty(:,sigt > 1e-12); 
% Vtesty = orth([ones(size(x(:))) Vtesty]); % ensure 1 in test space

disp('Done with Vtest svd')

rho = Vr'*U0(:,1);
rhou = Vr'*U0(:,2);
rhov = Vr'*U0(:,3);
E = Vr'*U0(:,4);

% surf(x,y,reshape(Vrp*rho,K,K))
% shading interp
% return

% tol = sqrt(sum(sum(dx*(U0 - Vr*[rho rhou rhov E]).^2)))/sqrt(sum(sum(dx*U0.^2)))

% hyperreduction
if 1                      
    
    wq0 = ones(size(x(:)))*dx^2;
    if 1
        
        Va = Vr;
        Vb = Vr;
        
%         [Vtest,Stest,~] = svd([Vtestx Vtesty],0);
%         Vb = Vtest(:,diag(Stest)>1e-10);
        
        sk = 1;
        Vmass = zeros(size(Va,1),size(Va,2)*size(Vb,2));
        for i = 1:size(Va,2)
            for j = 1:size(Vb,2)
                Vmass(:,sk) = Va(:,i).*Vb(:,j);
                sk = sk + 1;
            end
        end
        disp('Done building Vmass')
        
        % reduce target space        
        tic; [Vmass,Smass,~] = rsvd(Vmass,25*Nmodes); toc
        smass = diag(Smass);
        smass_energy = sqrt(1 - (cumsum(smass.^2)./sum(smass.^2)));
        Vmass = Vmass(:,smass_energy > tol);        
        disp('Done with Vmass svd')
                
        tic;[wr_unscaled id] = get_empirical_cubature(Vmass,ones(size(x(:))),tol,25*Nmodes); toc % scale to reduce roundoff
%         disp('Computing EQP points')
%         tic;[wr_unscaled id] = get_EQP_nodes(Vmass,ones(size(x(:))),tol,[],[]);toc
        wr = wr_unscaled*dx^2;
        
    else
        wr = wq0; id = 1:length(wq0);
    end    
    
    % check if test mass matrices are nonsingular
    Mtestx = Vtestx(id,:)'*diag(wr)*Vtestx(id,:);    
%     Ptestx = pinv(diag(sqrt(wr))*Vtestx(id,:))*diag(sqrt(wr));
%     Ptesty = pinv(diag(sqrt(wr))*Vtesty(id,:))*diag(sqrt(wr));
    
    if cond(Mtestx) > 1e8
        [Vx Dx] = eig(Mtestx);
        d = diag(Dx);
        [d,p] = sort(d,'ascend');
        Z = Vtestx*Vx(:,p(d < 1e-12)); 
        
        fprintf('initial cond of Mtest %g\n',cond(Mtestx))
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
        fprintf('\n getting additional points for stability\n')
        w0 = dx^2*ones(size(x(:)));
        [wz idz] = get_empirical_cubature(Zmass,w0,tol,size(Zmass,2));
        
        zscale = 1e-2;
        idt = unique([id; idz]);
        J = [Vmass zscale*Zmass];
        b = sum(J',2);
        w_unscaled = lsqlin(J(idt,:)',b,[],[],[] ,[] ,zeros(size(idt)),inf(size(idt)),J(idt,:)'\b); % upper bound = total vol
        w = w_unscaled*sum(w0)/K^2;
        
        Mtestx = Vtestx(idt,:)'*diag(w)*Vtestx(idt,:);
        %             [Vt Dt] = eig(Mtest);
        %             lam = diag(Dt);
        %             [minlam,idx] = min(lam);
        
        if cond(Mtestx) > 1e8
            keyboard
        end
        
        id = idt;
        wr = w;
        
        fprintf('Added %d new points, new cond of Mtest = %g\n',length(idz),cond(Mtestx))        
        Mtestx = Vtestx(id,:)'*diag(wr)*Vtestx(id,:);
    end        
    
    Mtesty = Vtesty(id,:)'*diag(wr)*Vtesty(id,:);    
    if cond(Mtesty) > 1e8
        keyboard
    end
    
    % make new mass, projection, interp matrices
    Mr = (Vr(id,:)'*diag(wr)*Vr(id,:));
    invM = inv(Mr);
    Pr = Mr\(Vr(id,:)'*diag(wr));
    Vr = Vr(id,:);
    
    Ptestx = Mtestx\(Vtestx(id,:)'*diag(wr));
    Ptesty = Mtesty\(Vtesty(id,:)'*diag(wr));
    
    Qx = Ptestx'*(Vtestx'*Qxfull*Vtestx)*Ptestx;
    Qy = Ptesty'*(Vtesty'*Qyfull*Vtesty)*Ptesty;    
    KS = Pr'*(Vrp'*Kfull*Vrp)*Pr;    
   
    KS = Vr'*KS;    
    VPr = Vr*Pr;
          
    % surface hyperreduction
    fid1 = find(abs(x-min(x(:)))<1e-10); 
    fid2 = find(abs(y-min(y(:)))<1e-10); 
    fid3 = find(abs(x-max(x(:)))<1e-10); 
    fid4 = find(abs(y-max(y(:)))<1e-10);
    fid = [fid1; fid2; fid3; fid4];    
    nx = [-ones(size(fid1)); zeros(size(fid2)); ones(size(fid3)); zeros(size(fid4))];
    ny = [zeros(size(fid2)); -ones(size(fid3)); zeros(size(fid4)); ones(size(fid1))];
    mapM = reshape(1:length(fid),K,4);
    mapP = [mapM(:,3); mapM(:,4); mapM(:,1); mapM(:,2)];        
       
    xf = x(fid);
    yf = y(fid);
    Vftestx = Vtestx(fid,:);  
    Vftesty = Vtesty(fid,:);  
    
%     for i = 1:length(mapP)
%         plot(xf,yf,'o')
%         hold on
%         plot(xf(i),yf(i),'x','linewidth',2,'markersize',14)
%         idP = mapP(i);
%         plot(xf(idP),yf(idP),'x','linewidth',2,'markersize',14)
%         pause
%     end
        
    Vf = Vrp(fid,:);
    VfPr = Vf*Pr;
       
    if 1
        Va = Vf;
        Vb = Vf;
        sk = 1;
        Vfmass = zeros(size(Va,1),size(Va,2)*size(Vb,2));
        for i = 1:size(Va,2)
            for j = 1:size(Vb,2)
                Vfmass(:,sk) = Va(:,i).*Vb(:,j);
                sk = sk + 1;
            end
        end
        disp('Done building Vfmass')
        
        [Vfmass,Smass,~] = rsvd(Vfmass,25*Nmodes);
        smass = diag(Smass);
        smass_energy = sqrt(1 - (cumsum(smass.^2)./sum(smass.^2)));
%         Vfmass = Vfmass(:,smass_energy > 1e-13); % in 1D
        Vfmass = Vfmass(:,smass_energy > tol); % in 1D
        disp('Done with Vfmass svd')
        
        % constrained lsq
        e = ones(size(Qxfull,2),1);
        C = [Vftestx'*diag(nx);
            Vftesty'*diag(ny)];
        d = [Vtestx'*Qxfull'*e; Vtesty'*Qyfull'*e];
        
        % linprog for constrained quadrature
        wf0 = dx*ones(size(Vftestx,1),1);
        
        disp('running linprog');       
        tic; [wf idf] = get_EQP_nodes(Vfmass,wf0,tol,C,d); toc        
        %tic; [wf_unscaled idf] = get_EQP_nodes(Vfmass,ones(size(Vftest,1),1),tol,C,d/dx); toc        
        %wf = wf_unscaled*dx;
    else
        wf = wf0; idf = 1:length(wf0);
    end
    
    norm(Vmass(id,:)'*wr - Vmass'*wq0)
    norm(Vfmass(idf,:)'*wf - Vfmass'*wf0)
    
    Efx = Vftestx(idf,:)*Ptestx;
    Efy = Vftesty(idf,:)*Ptesty;
    Bx = diag(nx(idf).*wf);
    By = diag(ny(idf).*wf);
    Z = zeros(length(wf));
    QNx = [Qx Efx'*Bx;
        -Bx*Efx Z];
%     BNx = [0*Qx 0*Efx'
%         0*Efx Bx];
    QNy = [Qy Efy'*By;
        -By*Efy Z];        
    
%     KS = Vr'*KS;
    
    Vf = Vf(idf,:);
    
    VN = [Vr;Vf];
    VfTW = Vf'*diag(wf);
    nx = nx(idf);
    ny = ny(idf);
    bxids = find(abs(abs(nx)-1) < 1e-10);
    byids = find(abs(abs(ny)-1) < 1e-10);
            
    clear Vtestx Vtesty Ptestx Ptesty

end

% return

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

res1 = zeros(Nmodes,1);
res2 = zeros(Nmodes,1);
res3 = zeros(Nmodes,1);
res4 = zeros(Nmodes,1);

rhs1 = zeros(Nmodes,1);
rhs2 = zeros(Nmodes,1);
rhs3 = zeros(Nmodes,1);
rhs4 = zeros(Nmodes,1);

% process blocks 
chunkSize = 32;
NT = size(QNx,1); % NT = total pts
Nchunk = ceil(NT/chunkSize);

for i = 1:Nsteps
    for INTRK = 1:5
        
        rhoq  = Vr*rho;
        rhouq = Vr*rhou;
        rhovq = Vr*rhov;
        Eq    = Vr*E;        
        
        v1 = Pr*V1(rhoq,rhouq,rhovq,Eq);
        v2 = Pr*V2(rhoq,rhouq,rhovq,Eq);
        v3 = Pr*V3(rhoq,rhouq,rhovq,Eq);
        v4 = Pr*V4(rhoq,rhouq,rhovq,Eq);
        
        v1N = Vr*v1;
        v2N = Vr*v2;
        v3N = Vr*v3;
        v4N = Vr*v4;

        rhoN  = U1(v1N,v2N,v3N,v4N);
        rhouN = U2(v1N,v2N,v3N,v4N);
        rhovN = U3(v1N,v2N,v3N,v4N);
        EN    = U4(v1N,v2N,v3N,v4N);
        uN    = rhouN./rhoN;
        vN    = rhovN./rhoN;
        
        % face vals 
        v1f = Vf*v1;
        v2f = Vf*v2;
        v3f = Vf*v3;
        v4f = Vf*v4;        
        rhoM  = U1(v1f,v2f,v3f,v4f);
        rhouM = U2(v1f,v2f,v3f,v4f);
        rhovM = U3(v1f,v2f,v3f,v4f);
        EM    = U4(v1f,v2f,v3f,v4f);
        uM    = rhouM./rhoM;
        vM    = rhovM./rhoM;
                
        % impose wall BCs
        rhoP = rhoM;
        uP = uM;
        vP = vM;
        EP = EM;
        uP(bxids) = -uM(bxids);
        vP(byids) = -vM(byids);
        
        % compute flux with LF penalty
        tau_lf = 1/2;
        unorm2 = (uM.^2+vM.^2);
        pM = (gamma-1)*(EM - .5*rhoM.*unorm2);
        cvel = sqrt(gamma*pM./rhoM);
        lam = sqrt(unorm2)+cvel;
        LFc = abs(lam);        
        % rhoP = rhoM, EP = EM so LF only activates 2nd penalty = mP= -mM
        LF1 = 0;
        LF2 = tau_lf*LFc.*(rhoM.*(-2*uM)).*wf; 
        LF3 = tau_lf*LFc.*(rhoM.*(-2*vM)).*wf; 
        LF4 = 0;
        
%         LF2 = 0;
%         LF3 = 0;
        
        f1 = nx.*fxS1(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + ny.*fyS1(rhoM,uM,vM,EM,rhoP,uP,vP,EP); 
        f2 = nx.*fxS2(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + ny.*fyS2(rhoM,uM,vM,EM,rhoP,uP,vP,EP) - LF2;  
        f3 = nx.*fxS3(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + ny.*fyS3(rhoM,uM,vM,EM,rhoP,uP,vP,EP) - LF3;  
        f4 = nx.*fxS4(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + ny.*fyS4(rhoM,uM,vM,EM,rhoP,uP,vP,EP);  
                    
        % stack vol/face vars
        rhoN = [rhoN;rhoM];
        uN = [uN;uM];
        vN = [vN;vM];
        EN = [EN;EM];
        betaN = beta(rhoN,uN,vN,EN);
        
        r1 = zeros(NT,1);
        r2 = zeros(NT,1);
        r3 = zeros(NT,1);
        r4 = zeros(NT,1);
        for ii=1:Nchunk
            for jj = 1:Nchunk
                id1 = (1:chunkSize) + (ii-1)*chunkSize;
                id2 = (1:chunkSize) + (jj-1)*chunkSize;
                if max(id1)> NT
                    sz = NT-(ii-1)*chunkSize;
                    id1 = (1:sz) + (ii-1)*chunkSize;
                end
                if max(id2) > NT
                    sz = NT-(jj-1)*chunkSize;
                    id2 = (1:sz) + (jj-1)*chunkSize;
                end
                
                Qxij = QNx(id1,id2);
                Qyij = QNy(id1,id2);
                
                [rhox rhoy] = meshgrid(rhoN(id1),rhoN(id2));
                [ux uy]     = meshgrid(uN(id1),uN(id2));
                [vx vy]     = meshgrid(vN(id1),vN(id2));
                %[Ex Ey]     = meshgrid(EN(id1),EN(id2));
                [betax betay]     = meshgrid(betaN(id1),betaN(id2));                                
                
                [FxS1 FxS2 FxS3 FxS4 FyS1 FyS2 FyS3 FyS4] = ...
                    euler_flux_2d(rhox,rhoy,ux,uy,vx,vy,betax,betay);
                r1(id1) = r1(id1) + (sum(Qxij.*FxS1',2) + sum(Qyij.*FyS1',2));
                r2(id1) = r2(id1) + (sum(Qxij.*FxS2',2) + sum(Qyij.*FyS2',2));
                r3(id1) = r3(id1) + (sum(Qxij.*FxS3',2) + sum(Qyij.*FyS3',2));
                r4(id1) = r4(id1) + (sum(Qxij.*FxS4',2) + sum(Qyij.*FyS4',2));

%                 r1(id1) = r1(id1) + (sum(Qxij.*fxS1(rhox,ux,vx,Ex,rhoy,uy,vy,Ey)',2) + sum(Qyij.*fyS1(rhox,ux,vx,Ex,rhoy,uy,vy,Ey)',2));
%                 r2(id1) = r2(id1) + (sum(Qxij.*fxS2(rhox,ux,vx,Ex,rhoy,uy,vy,Ey)',2) + sum(Qyij.*fyS2(rhox,ux,vx,Ex,rhoy,uy,vy,Ey)',2));
%                 r3(id1) = r3(id1) + (sum(Qxij.*fxS3(rhox,ux,vx,Ex,rhoy,uy,vy,Ey)',2) + sum(Qyij.*fyS3(rhox,ux,vx,Ex,rhoy,uy,vy,Ey)',2));
%                 r4(id1) = r4(id1) + (sum(Qxij.*fxS4(rhox,ux,vx,Ex,rhoy,uy,vy,Ey)',2) + sum(Qyij.*fyS4(rhox,ux,vx,Ex,rhoy,uy,vy,Ey)',2));
            end
        end                

        r1 = VN'*r1  + VfTW * f1;
        r2 = VN'*r2  + VfTW * f2;
        r3 = VN'*r3  + VfTW * f3;
        r4 = VN'*r4  + VfTW * f4;
        
        % convective stability test
        rhstest(i) = sum(sum(v1.*r1 + v2.*r2 + v3.*r3 + v4.*r4));                
        
        % add dissipation + invert mass                
        rhs1 = -invM*(r1 + tau*KS*rhoq);
        rhs2 = -invM*(r2 + tau*KS*rhouq);
        rhs3 = -invM*(r3 + tau*KS*rhovq);
        rhs4 = -invM*(r4 + tau*KS*Eq   );
        
%         keyboard
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        res4 = rk4a(INTRK)*res4 + dt*rhs4;
        
        rho  = rho   + rk4b(INTRK)*res1;
        rhou = rhou  + rk4b(INTRK)*res2;
        rhov = rhov  + rk4b(INTRK)*res3;
        E    = E     + rk4b(INTRK)*res4;
    end
    
    Uhi = reshape(Vrp*[rho rhou rhov E],K*K,4); 
    rhoROM = Uhi(:,1);
    rhouROM = Uhi(:,2);
    rhovROM = Uhi(:,3);
    EROM = Uhi(:,4);    
    entropyROM(i) = sum(sum(-rhoROM.*sfun(rhoROM,rhouROM,rhovROM,EROM)))*dx^2;
    
    Uhi = reshape(Usnap(:,i),K*K,4); 
    rhoFOM = Uhi(:,1);
    rhouFOM = Uhi(:,2);
    rhovFOM = Uhi(:,3);
    EFOM = Uhi(:,4);
    entropyFOM(i) = sum(sum(-rhoFOM.*sfun(rhoFOM,rhouFOM,rhovFOM,EFOM)))*dx^2;
    
    if mod(i,5) == 0 || i==Nsteps
        surf(x,y,reshape(Vrp*rho,K,K))
        shading interp
        view(2)
        
        Usnapi = reshape(Usnap(:,i),K*K,4);
        Uhi = Vrp*[rho rhou rhov E];
        e = Uhi - Usnapi;
        err = sqrt(sum(sum(dx^2*(e.^2))))/sqrt(sum(sum(dx^2*Uhi.^2)));
        title(sprintf('time = %f, step %d / %d, rhstest = %g, err = %g\n',dt*i,i,Nsteps,rhstest(i),err));
        drawnow
    end
    
end

figure
minS = min(min(entropyROM),min(entropyFOM));
semilogy(entropyFOM-minS,'-','linewidth',2)
hold on
semilogy(entropyROM-minS,'--','linewidth',2)

%%
Usnapi = reshape(Usnap(:,i),K*K,4);
Uhi = Vrp*[rho rhou rhov E];
e = Uhi - Usnapi;

figure
surf(x,y,reshape(Usnapi(:,1),K,K)); 
shading interp; view(2)
axis off; axis equal
% print(gcf,'-dpng','~/Desktop/bbwadg/docs/ESDG_ROM/figs/pulse2d.png')

figure
surf(x,y,reshape(Uhi(:,1),K,K));
shading interp; view(2)
axis off; axis equal; 

% print(gcf,'-dpng','~/Desktop/bbwadg/docs/ESDG_ROM/figs/pulse2d_ROM.png')
figure
surf(x,y,reshape(Uhi(:,1),K,K));
shading interp; view(2)
axis off; axis equal; 

hold on
plot3(x(id),y(id),Vr*rho+.1,'r.','linewidth',2,'markersize',12)
plot3(xf(idf),yf(idf),Vf*rho+.1,'b.','linewidth',2,'markersize',12)

% print(gcf,'-dpng','~/Desktop/bbwadg/docs/ESDG_ROM/figs/pulse2d_ROM_pts.png')

% err = sqrt(sum(sum(dx^2*(e.^2))))/sqrt(sum(sum(dx^2*Uhi.^2)));
% title(sprintf('err = %g\n',err))
