clear
K = 200;
%FinalTime = .35;
FinalTime = 2.5;

xv = linspace(-1,1,K+1)';
x = .5*(xv(1:end-1)+xv(2:end));
dx = (max(xv(:))-min(xv(:)))/K;

tau = .1*dx;

e = ones(K-1,1);
S = diag(e,1)-diag(e,-1);
S(1,:) = 0; S(:,1) = 0;
S(end,:) = 0; S(:,end) = 0;
S(1,end) = -1; S(end,1) = 1;
S(1,2) = 1; S(2,1) = -1;
S(K-1,K) = 1; S(K,K-1) = -1;
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

rho = ones(size(x)) + exp(-50*((x+.5).^2+(y+.5).^2));
% rho = 1 + (abs(x)<.5).*(abs(y)<.5);
% xc = x + .5; yc = y + .5;
% rho = 2*ones(size(x)) + .5*exp(-100*(xc.^2+0*yc.^2)); %.*sin(pi*xc).*sin(pi*yc);
u = zeros(size(x));
v = zeros(size(x));
p = rho.^(gamma);

rhou = rho.*u;
rhov = rho.*v;
E    = p/(gamma-1) + .5*rho.*(u.^2+v.^2);

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

% KH instab
a = .1; sig = .1; %5*sqrt(2)*1e-3;
ff = 1./(1+exp(-(y+.5)./(sig.^2)))-1./(1+exp(-(y-.5)./(sig.^2)));
rho = ff + 1;
v = a*sin(2*pi*x).*(exp(-(y-.5).^2/((sig/2).^2))+exp(-(y+.5).^2/((sig/2).^2)));
u = ff - .5;
p = 2.5*ones(size(x));

rhou = rho.*u;
rhov = rho.*v;
E    = p/(gamma-1) + .5*rho.*(u.^2+v.^2);
% surf(x,y,rho);return


%% make ROM ops

load Usnap_euler2d_kh.mat

U0 = reshape(Usnap(:,1),K^2,4);
Us1 = Usnap(1:K^2,:);
Us2 = Usnap((1:K^2) + K^2,:);
Us3 = Usnap((1:K^2) + 2*K^2,:);
Us4 = Usnap((1:K^2) + 3*K^2,:);

[Vr, Sr, ~] = svd([Us1 Us2 Us3 Us4 V1(Us1,Us2,Us3,Us4) V2(Us1,Us2,Us3,Us4) V3(Us1,Us2,Us3,Us4) V4(Us1,Us2,Us3,Us4)],0);
sig = diag(Sr);
disp('Done with Vr svd')
% semilogy(sig,'o')
keyboard

Nmodes = 75;
tol = sqrt(sum(sig(Nmodes+1:end).^2)/sum(sig.^2))
% tol = sqrt(sum(sum(dx*(U0 - Vr*[rho rhou rhov E]).^2)))/sqrt(sum(sum(dx*U0.^2)))

Vr = Vr(:,1:Nmodes);
Vrp = Vr;

% add range + constants to snapshots
Qxfull = kron(S,dx*speye(K));
Qyfull = kron(dx*speye(K),S);
Kfull = kron(KS,dx*speye(K)) + kron(dx*speye(K),KS); % d2u/dx2 + d2u/dy2;
%DVr = [Qxfull*Vr Qyfull*Vr];
DVrx = [Qxfull*Vr];
DVry = [Qyfull*Vr];
% DVr = [Qyfull*Vr]; %TESTING ONLY


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

[Vtestx, Stest, ~] = svd([Vr DVrx],0);toc
sigt = diag(Stest);
Vtestx = Vtestx(:,sigt > 1e-12); 
Vtestx = orth([ones(size(x(:))) Vtestx]); % ensure 1 in test space

[Vtesty, Stest, ~] = svd([Vr DVry],0);toc
sigt = diag(Stest);
Vtesty = Vtesty(:,sigt > 1e-12); 
Vtesty = orth([ones(size(x(:))) Vtesty]); % ensure 1 in test space

disp('Done with Vtest svd')

% hyperreduction
if 1
    
    Va = Vr;
    Vb = Vr; 
%     Vb = Vtest;
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
    
    % reduce target space
%     [Vmass,Smass,~] = svd(Vmass,0);    
    tic; [Vmass,Smass,~] = rsvd(Vmass,1500); toc
    smass = diag(Smass);
    smass_energy = sqrt(1 - (cumsum(smass.^2)./sum(smass.^2)));
    Vmass = Vmass(:,smass_energy > tol);
%     Vmass = Vmass(:,1:4*Nmodes);    
    disp('Done with Vmass svd')
        
%     [wr id] = get_empirical_cubature(Vmass,ones(size(x(:)))*dx^2,tol,25*Nmodes); % cheaper
    [wr_unscaled id] = get_empirical_cubature(Vmass,ones(size(x(:))),tol,25*Nmodes); % don't scale to reduce roundoff
    wr = wr_unscaled*dx^2;        
    
    % make new mass, projection, interp matrices
    Mr = (Vrp(id,:)'*diag(wr)*Vrp(id,:));
    invM = inv(Mr);
    Pr = Mr\(Vrp(id,:)'*diag(wr));
    Vr = Vrp(id,:);    
    
%     Mtest = Vtest(id,:)'*diag(wr)*Vtest(id,:);
%     Ptest = Mtest\(Vtest(id,:)'*diag(wr));
    
    Mtestx = Vtestx(id,:)'*diag(wr)*Vtestx(id,:);
    Mtesty = Vtesty(id,:)'*diag(wr)*Vtesty(id,:);
    
    if cond(Mtestx) > 1e8
        keyboard
    end
    
    if cond(Mtesty) > 1e8
        keyboard
    end
    
    Ptestx = Mtestx\(Vtestx(id,:)'*diag(wr));
    Ptesty = Mtesty\(Vtesty(id,:)'*diag(wr));
    
    Qx = Ptestx'*(Vtestx'*Qxfull*Vtestx)*Ptestx;
    Qy = Ptesty'*(Vtesty'*Qyfull*Vtesty)*Ptesty;    
    KS = Pr'*(Vrp'*Kfull*Vrp)*Pr;    
   
    KS = Vr'*KS;    
    VPr = Vr*Pr;
    
    clear Vtestx Vtesty Ptestx Ptesty
    
%     keyboard
%     if max(dx*dx./eig(Mr)) > 2
%         % if CFL = bad
%         keyboard
%     end
%     
%     keyboard

    
end

rho  = Vrp'*U0(:,1);
rhou = Vrp'*U0(:,2);
rhov = Vrp'*U0(:,3);
E    = Vrp'*U0(:,4);

% surf(x,y,reshape(Vrp*rho,K,K))
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


dt = 1*dx;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

res1 = zeros(size(rho));
res2 = zeros(size(rho));
res3 = zeros(size(rho));
res4 = zeros(size(rho));

rhs1 = zeros(size(rho));
rhs2 = zeros(size(rho));
rhs3 = zeros(size(rho));
rhs4 = zeros(size(rho));

chunkSize = 32;
Nchunk = ceil(size(Vr,1)/chunkSize);

for i = 1:Nsteps
    for INTRK = 1:5       
        
        rhoq  = Vr*rho;
        rhouq = Vr*rhou;
        rhovq = Vr*rhov;
        Eq    = Vr*E;
        
        v1N = VPr*V1(rhoq,rhouq,rhovq,Eq);
        v2N = VPr*V2(rhoq,rhouq,rhovq,Eq);
        v3N = VPr*V3(rhoq,rhouq,rhovq,Eq);
        v4N = VPr*V4(rhoq,rhouq,rhovq,Eq);
        
        rhoN  = U1(v1N,v2N,v3N,v4N);
        rhouN = U2(v1N,v2N,v3N,v4N);
        rhovN = U3(v1N,v2N,v3N,v4N);
        EN    = U4(v1N,v2N,v3N,v4N);
        uN    = rhouN./rhoN;
        vN    = rhovN./rhoN;
        
        betaN = beta(rhoN,uN,vN,EN);
        
        r1 = zeros(size(Vr,1),1);
        r2 = zeros(size(Vr,1),1);
        r3 = zeros(size(Vr,1),1);
        r4 = zeros(size(Vr,1),1);
        for ii=1:Nchunk
            for jj = 1:Nchunk
                id1 = (1:chunkSize) + (ii-1)*chunkSize;
                id2 = (1:chunkSize) + (jj-1)*chunkSize;
                if max(id1)>size(Vr,1)
                    sz = size(Vr,1)-(ii-1)*chunkSize;
                    id1 = (1:sz) + (ii-1)*chunkSize;
                end
                if max(id2)>size(Vr,1)
                    sz = size(Vr,1)-(jj-1)*chunkSize;
                    id2 = (1:sz) + (jj-1)*chunkSize;                    
                end
                
                Qxij = Qx(id1,id2);
                Qyij = Qy(id1,id2);
                
                [rhox rhoy] = meshgrid(rhoN(id1),rhoN(id2));
                [ux uy]     = meshgrid(uN(id1),uN(id2));
                [vx vy]     = meshgrid(vN(id1),vN(id2));
                %[Ex Ey]     = meshgrid(EN(id1),EN(id2));
                [betax betay]     = meshgrid(betaN(id1),betaN(id2));
                
                [FxS1 FxS2 FxS3 FxS4 FyS1 FyS2 FyS3 FyS4] = euler_flux_2d(rhox,rhoy,ux,uy,vx,vy,betax,betay)
                r1(id1) = r1(id1) + (sum(Qxij.*FxS1',2) + sum(Qyij.*fyS1',2));
                r2(id1) = r2(id1) + (sum(Qxij.*fxS2',2) + sum(Qyij.*fyS2',2));
                r3(id1) = r3(id1) + (sum(Qxij.*fxS3',2) + sum(Qyij.*fyS3',2));
                r4(id1) = r4(id1) + (sum(Qxij.*fxS4',2) + sum(Qyij.*fyS4',2));
                
%                 r1(id1) = r1(id1) + (sum(Qxij.*fxS1(rhox,ux,vx,Ex,rhoy,uy,vy,Ey)',2) + sum(Qyij.*fyS1(rhox,ux,vx,Ex,rhoy,uy,vy,Ey)',2));
%                 r2(id1) = r2(id1) + (sum(Qxij.*fxS2(rhox,ux,vx,Ex,rhoy,uy,vy,Ey)',2) + sum(Qyij.*fyS2(rhox,ux,vx,Ex,rhoy,uy,vy,Ey)',2));
%                 r3(id1) = r3(id1) + (sum(Qxij.*fxS3(rhox,ux,vx,Ex,rhoy,uy,vy,Ey)',2) + sum(Qyij.*fyS3(rhox,ux,vx,Ex,rhoy,uy,vy,Ey)',2));
%                 r4(id1) = r4(id1) + (sum(Qxij.*fxS4(rhox,ux,vx,Ex,rhoy,uy,vy,Ey)',2) + sum(Qyij.*fyS4(rhox,ux,vx,Ex,rhoy,uy,vy,Ey)',2));
            end
        end
        
        % convective stability test
        rhstest(i) = sum(sum(v1N.*r1 + v2N.*r2 + v3N.*r3 + v4N.*r4));
        
        % add dissipation + invert mass
        rhs1 = -invM*(Vr'*r1 + tau*KS*rhoq);
        rhs2 = -invM*(Vr'*r2 + tau*KS*rhouq);
        rhs3 = -invM*(Vr'*r3 + tau*KS*rhovq);
        rhs4 = -invM*(Vr'*r4 + tau*KS*Eq);
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        res4 = rk4a(INTRK)*res4 + dt*rhs4;
        
        rho  = rho   + rk4b(INTRK)*res1;
        rhou = rhou  + rk4b(INTRK)*res2;
        rhov = rhov  + rk4b(INTRK)*res3;
        E    = E     + rk4b(INTRK)*res4;
        
    end
    
    if mod(i,5) == 0
        surf(x,y,reshape(Vrp*rho,K,K))
        shading interp
        view(2)
        title(sprintf('time = %f, step %d / %d, rhstest = %g\n',dt*i,i,Nsteps,rhstest(i)));
        drawnow
    end
end

sqrt(sum(sum(dx^2*(Vrp*rho-Us1(:,end)).^2)))/sqrt(sum(sum(dx^2*Us1(:,end))))
