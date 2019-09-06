clear
K = 25;
FinalTime = .3;
xv = linspace(-1,1,K+1)';
x1D = .5*(xv(1:end-1)+xv(2:end));
x = x1D;
dx = 1/K;

e = ones(K-1,1);
S = diag(e,1)-diag(e,-1);
S(1,:) = 0; S(:,1) = 0;
S(end,:) = 0; S(:,end) = 0;
% S(1,end) = -1; S(end,1) = 1;
S(1,2) = 1; S(2,1) = -1;
S(K-1,K) = 1; S(K,K-1) = -1;

Q = S;
Q(1) = -1;
Q(end) = 1; 
% Q = .5*Q;
% M = dx/2 * eye(K);

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

%% make ROM

load Usnap_euler2d_wall

U0 = reshape(Usnap(:,1),K^2,4);
u1 = Usnap(1:K^2,:);
u2 = Usnap((1:K^2) + K^2,:);
u3 = Usnap((1:K^2) + 2*K^2,:);
u4 = Usnap((1:K^2) + 3*K^2,:);

[Vr, Sr, ~] = svd([u1 u2 u3 u4 ...
    V1(u1,u2,u3,u4) V2(u1,u2,u3,u4) V3(u1,u2,u3,u4) V4(u1,u2,u3,u4)],0); 
sig = diag(Sr); 
% semilogy(sig,'o')

Nmodes = 25;

Vr = Vr(:,1:Nmodes);
Vrp = Vr;

tol = sum(sig(Nmodes+1:end).^2)/sum(sig.^2)

% add range + constants to snapshots
Qxfull = kron(S,dx*speye(K));
Qyfull = kron(dx*speye(K),S);
Kfull = kron(KS,dx*speye(K)) + kron(dx*speye(K),KS); % d2u/dx2 + d2u/dy2;

[Vtest, Stest, ~] = svd([Vr Qxfull*Vr Qyfull*Vr]);
sigt = diag(Stest);
Vtest = Vtest(:,sigt > 1e-12);
Vtest = orth([ones(size(x(:))) Vtest]); % ensure 1 in test space

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
           
    Va = Vr;
%     Vb = Vr; 
    Vb = Vtest;
    sk = 1;
    Vmass = zeros(size(Va,1),size(Va,2)*size(Vb,2));
    for i = 1:size(Va,2)
        for j = 1:size(Vb,2)
            Vmass(:,sk) = Va(:,i).*Vb(:,j);
            sk = sk + 1;
        end
    end
    disp('Done building Vmass')
%     keyboard
    
%     % reduce target space
%     tic; [Vmass,Smass,~] = rsvd(Vmass,25*Nmodes); toc
%     smass = diag(Smass);
%     smass_energy = sqrt(1 - (cumsum(smass.^2)./sum(smass.^2)));
%     Vmass = Vmass(:,smass_energy > tol);
% %     Vmass = Vmass(:,1:4*Nmodes);    
%     disp('Done with Vmass svd')
        
    wq0 = ones(size(x(:)))*dx^2;    

%     [wr id] = get_empirical_cubature(Vmass,ones(size(x(:)))*dx^2,tol,25*Nmodes); % cheaper
%     [wr_unscaled id] = get_empirical_cubature(Vmass,ones(size(x(:))),tol,25*Nmodes); % don't scale to reduce roundoff
%     wr = wr_unscaled*dx^2;        
%     tic;[wr id] = get_EQP_nodes(Vmass,wq0,tol,[],[]);toc    

   
    wr = wq0; id = 1:length(wq0);

    % make new mass, projection, interp matrices
    Mr = (Vr(id,:)'*diag(wr)*Vr(id,:));
    invM = inv(Mr);
    Pr = Mr\(Vr(id,:)'*diag(wr));
    Vr = Vr(id,:);
    
    Mtest = Vtest(id,:)'*diag(wr)*Vtest(id,:);
    Ptest = Mtest\(Vtest(id,:)'*diag(wr));
    Qx = Ptest'*(Vtest'*Qxfull*Vtest)*Ptest;
    Qy = Ptest'*(Vtest'*Qyfull*Vtest)*Ptest;
    KS = Ptest'*(Vtest'*Kfull*Vtest)*Ptest; 
    
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
    Vftest = Vtest(fid,:);  
    
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
       
    Va = Vftest;
    Vb = Vftest; 
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
    Vfmass = Vfmass(:,smass_energy > 1e-13); % in 1D
    disp('Done with Vfmass svd')
    
    % constrained lsq
    e = ones(size(Qxfull,2),1);
    C = [Vftest'*diag(nx);
        Vftest'*diag(ny)];
    d = [Vtest'*Qxfull'*e; Vtest'*Qyfull'*e];
    
    % linprog for constrained quadrature
    wf0 = dx*ones(size(Vftest,1),1);            
    
    f = ones(size(wf0));
        
%     disp('running linprog')
%     tic;
%     [wf idf] = get_EQP_nodes(Vfmass,wf0,tol,C,d);
%     toc
    wf = wf0; idf = 1:length(wf0);
    
    norm(Vmass(id,:)'*wr - Vmass'*wq0)
    norm(Vfmass(idf,:)'*wf - Vfmass'*wf0)
    
    Ef = Vftest(idf,:)*Ptest;
    Bx = diag(nx(idf).*wf);
    By = diag(ny(idf).*wf);
    Z = zeros(length(wf));
    QNx = [Qx Ef'*Bx;
        -Bx*Ef Z];
    BNx = [0*Qx 0*Ef'
        0*Ef Bx];
    BNy = [0*Qx 0*Ef'
        0*Ef By];
    QNy = [Qy Ef'*By;
        -By*Ef Z];        
    VN = [Vr;Vf];
    
    KS = Vr'*KS;
    
    Vf = Vf(idf,:);
    VfTW = Vf'*diag(wf);
    nx = nx(idf);
    ny = ny(idf);
    bxids = find(abs(abs(nx)-1) < 1e-10);
    byids = find(abs(abs(ny)-1) < 1e-10);
            
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
         
     
dt = .5*dx;
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
        uM    = U2(v1f,v2f,v3f,v4f)./rhoM;
        vM    = U3(v1f,v2f,v3f,v4f)./rhoM;
        EM    = U4(v1f,v2f,v3f,v4f);        
        
        % impose wall BCs
        rhoP = rhoM(mapP);
        uP = uM(mapP);
        vP = vM(mapP);
        EP = EM(mapP);        
        
%         % impose wall BCs
%         rhoP = rhof;
%         uP = uf;
%         vP = vf;
%         EP = Ef;
%         uP(bxids) = -uf(bxids);
%         vP(byids) = -vf(byids);
        
        % compute flux
        % add LF penalty
        f1 = nx.*fxS1(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + ny.*fyS1(rhoM,uM,vM,EM,rhoP,uP,vP,EP); 
        f2 = nx.*fxS2(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + ny.*fyS2(rhoM,uM,vM,EM,rhoP,uP,vP,EP);  
        f3 = nx.*fxS3(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + ny.*fyS3(rhoM,uM,vM,EM,rhoP,uP,vP,EP);  
        f4 = nx.*fxS4(rhoM,uM,vM,EM,rhoP,uP,vP,EP) + ny.*fyS4(rhoM,uM,vM,EM,rhoP,uP,vP,EP);  
            
        tau = 0*.5*dx;
        
        % stack vol/face vars
        rhoN = [rhoN;rhoM];
        uN = [uN;uM];
        vN = [vN;vM];
        EN = [EN;EM];               
        
        [rhox rhoy] = meshgrid(rhoN);
        [ux uy]     = meshgrid(uN);
        [vx vy]     = meshgrid(vN);
        [Ex Ey]     = meshgrid(EN);
        
        Fx1 = fxS1(rhox,ux,vx,Ex,rhoy,uy,vy,Ey);
        Fx2 = fxS2(rhox,ux,vx,Ex,rhoy,uy,vy,Ey);
        Fx3 = fxS3(rhox,ux,vx,Ex,rhoy,uy,vy,Ey);
        Fx4 = fxS4(rhox,ux,vx,Ex,rhoy,uy,vy,Ey);
        Fy1 = fyS1(rhox,ux,vx,Ex,rhoy,uy,vy,Ey);
        Fy2 = fyS2(rhox,ux,vx,Ex,rhoy,uy,vy,Ey);
        Fy3 = fyS3(rhox,ux,vx,Ex,rhoy,uy,vy,Ey);
        Fy4 = fyS4(rhox,ux,vx,Ex,rhoy,uy,vy,Ey);
        
        r1 = VN'*(sum(QNx.*Fx1,2)+sum(QNy.*Fy1,2)) + VfTW * f1;
        r2 = VN'*(sum(QNx.*Fx2,2)+sum(QNy.*Fy2,2)) + VfTW * f2;
        r3 = VN'*(sum(QNx.*Fx3,2)+sum(QNy.*Fy3,2)) + VfTW * f3;
        r4 = VN'*(sum(QNx.*Fx4,2)+sum(QNy.*Fy4,2)) + VfTW * f4;
        
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
    
    if mod(i,1) == 0
        surf(x,y,reshape(Vrp*rho,K,K))
        shading interp
        view(2)
        title(sprintf('time = %f, step %d / %d, rhstest = %g\n',dt*i,i,Nsteps,rhstest(i)));
        drawnow
    end
end
