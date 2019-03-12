% SBP spectral method

clear

N = 100; % must be even
FinalTime = 4;

Np1D = N;

dx = 2/N;
x = linspace(-1,1,N+1); x = x(1:end-1)';

h = 2*pi/N; % original h on [-pi,pi]
column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)];
D = toeplitz(column,column([1 N:-1:2]))*pi;

W = dx*eye(Np1D);
Q = W*D;
K = D'*W*D;


%% make 2D
[x y] = meshgrid(x);
Np = Np1D^2;

%%

% Euler
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

%%

p = ones(size(x));
width = .1;
sig = 10000;
ff = 1-(1./(1+exp(-sig*(y-width))) + 1./(1+exp(sig*(y+width))));
dv = .1*sin(2*pi*x).*ff;
du = .5*ff;
v = zeros(size(x)) + dv;
u = -.5*ones(size(x)) + du;

rho = ones(size(x));
rhou = rho.*u;
rhov = rho.*v;
E    = p/(gamma-1) + .5*rho.*(u.^2+v.^2);

% surf(x,y,rhov);shading interp;view(2);return

%% build reduced model

load Usnap_euler_spectral

U0 = reshape(Usnap(:,1),Np,4);
u1 = Usnap((1:Np),:);
u2 = Usnap((1:Np)+Np,:);
u3 = Usnap((1:Np)+2*Np,:);
u4 = Usnap((1:Np)+3*Np,:);
Us = [u1 u2 u3 u4 V1(u1,u2,u3,u4) V2(u1,u2,u3,u4) V3(u1,u2,u3,u4) V4(u1,u2,u3,u4)];
[Vr, Sr, ~] = svd(Us,0);
sig = diag(Sr);

Nmodes = 35;
Vr = Vr(:,1:Nmodes);
Vrp = Vr;

% est tol based on energy
tol = sqrt(sum(sig(Nmodes+1:end).^2)./sum(sig.^2));

% full ops
Qxfull = kron(Q,W);
Qyfull = kron(W,Q);
Kfull = kron(K,W) + kron(W,K); % d2u/dx2 + d2u/dy2

% test space
[Vtest, Stest, ~] = svd([Vr Qxfull*Vr Qyfull*Vr],0);
sigt = diag(Stest);
sigerr = sqrt(1-cumsum(sigt.^2)/sum(sigt.^2));
Vtest = Vtest(:,sigerr > 1e-8);
% Vtest = Vtest(:,1:25);

% ensure constant in space
Vtest = orth([ones(Np,1) Vtest]); 

rho = Vr'*U0(:,1);
rhou = Vr'*U0(:,2);
rhov = Vr'*U0(:,3);
E = Vr'*U0(:,4);

err0 = sqrt(sum(dx^2*(Vrp*rhou-U0(:,2)).^2))/sqrt(sum(dx^2*(Vrp*rhou).^2));
fprintf('tol = %g, initial err = %g\n',tol,err0)

% hyperreduc
if 1
%     Vtarget = Vtest;
%     Vmass = zeros(Np,size(Vtarget,2)*(size(Vtarget,2)+1)/2);
%     sk = 1;
%     for i = 1:size(Vtarget,2)
%         for j = i:size(Vtarget,2)
%             Vmass(:,sk) = Vtarget(:,i).*Vtarget(:,j);
%             sk = sk + 1;
%         end
%     end
%     b = sum(Vmass'*dx^2,2);

    Vmass = zeros(Np,size(Vr,2)*size(Vtest,2));
    sk = 1;
    for i = 1:size(Vr,2)
        for j = i:size(Vtest,2)
            Vmass(:,sk) = Vr(:,i).*Vtest(:,j);
            sk = sk + 1;
        end
    end
    b = sum(Vmass'*dx^2,2);    
    
    maxpts = Nmodes*10;
%     id = get_empirical_cubature(Vmass,b,tol,maxpts);
    id = get_empirical_cubature(Vmass,b,err0,inf);
%     id = get_empirical_cubature(Vmass,b,tol,inf);
    wr = Vmass(id,:)'\b;    
        
    % make new ops
    Mtest = Vtest(id,:)'*diag(wr)*Vtest(id,:);
    Ptest = Mtest\(Vtest(id,:)'*diag(wr));    
    Qx = Ptest'*(Vtest'*Qxfull*Vtest)*Ptest;
    Qy = Ptest'*(Vtest'*Qyfull*Vtest)*Ptest;
    K = Ptest'*(Vtest'*Kfull*Vtest)*Ptest;
    
    % make new mass, projection, interp matrices
    Mr = (Vr(id,:)'*diag(wr)*Vr(id,:));    
    Pr = Mr\(Vr(id,:)'*diag(wr));
    invMVrT = Mr\(Vr(id,:)');
    Vr = Vr(id,:);
    
else
    
    % full ops
    Qx = Vtest*(Vtest'*Qxfull*Vtest)*Vtest';
    Qy = Vtest*(Vtest'*Qyfull*Vtest)*Vtest';
    K = Vtest*(Vtest'*Kfull*Vtest)*Vtest';
    invM = Vr'*Vr*(1/dx^2);
    Pr = invM*(Vr'*dx^2);
    
end


%%

dt = dx;
% dt = 2*dt;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;


u = 0*rho;
res1 = zeros(size(u));
res2 = zeros(size(u));
res3 = zeros(size(u));
res4 = zeros(size(u));
rhs1 = zeros(size(u));
rhs2 = zeros(size(u));
rhs3 = zeros(size(u));
rhs4 = zeros(size(u));

for i = 1:Nsteps
    for INTRK = 1:5                

        rhoq = Vr*rho;
        rhouq = Vr*rhou;
        rhovq = Vr*rhov;
        Eq = Vr*E;
        
        % project entropy vars
        v1 = Pr*V1(rhoq,rhouq,rhovq,Eq);
        v2 = Pr*V2(rhoq,rhouq,rhovq,Eq);
        v3 = Pr*V3(rhoq,rhouq,rhovq,Eq);
        v4 = Pr*V4(rhoq,rhouq,rhovq,Eq);        
        v1q = Vr*v1;
        v2q = Vr*v2;
        v3q = Vr*v3;
        v4q = Vr*v4;
        
        % entropy projection
        rhoq = U1(v1q,v2q,v3q,v4q);
        rhouq = U2(v1q,v2q,v3q,v4q);
        rhovq = U3(v1q,v2q,v3q,v4q);
        Eq = U4(v1q,v2q,v3q,v4q);
        
        % eval interactions
        uq = rhouq./rhoq;
        vq = rhovq./rhoq;
        betaq = beta(rhoq,uq,vq,Eq);        
        [rhox rhoy] = meshgrid(rhoq);
        [ux uy] = meshgrid(uq);
        [vx vy] = meshgrid(vq);
        [betax betay] = meshgrid(betaq);
        
        % aux terms
        rholog = logmean(rhox,rhoy);
        rhoavg = avg(rhox,rhoy);
        uavg = avg(ux,uy);
        vavg = avg(vx,vy);
        vnavg = 2*(uavg.^2 + vavg.^2) - (avg(ux.^2,uy.^2) + avg(vx.^2,vy.^2));
        pa = rhoavg./(2*avg(betax,betay));
        
        % x fluxes
        FxS1 = rholog.*uavg;
        FxS2 = FxS1.*uavg + pa;
        FxS3 = FxS1.*vavg;
        f4aux = rholog./(2*(gamma-1)*logmean(betax,betay)) + pa + .5*rholog.*vnavg;
        FxS4 = f4aux.*uavg;
        
        % y fluxes
        FyS1 = rholog.*vavg;
        FyS2 = FyS1.*uavg;
        FyS3 = FyS1.*vavg + pa;
        f4aux = rholog./(2*(gamma-1)*logmean(betax,betay)) + pa + .5*rholog.*vnavg;
        FyS4 = f4aux.*vavg;
        
        % flux differencing
        rhs1 = sum(Qx.*FxS1,2) + sum(Qy.*FyS1,2);
        rhs2 = sum(Qx.*FxS2,2) + sum(Qy.*FyS2,2);
        rhs3 = sum(Qx.*FxS3,2) + sum(Qy.*FyS3,2);
        rhs4 = sum(Qx.*FxS4,2) + sum(Qy.*FyS4,2);        
        
        % add diffusion and project down
        tau = 5e-4;
        rhs1 = -invMVrT*(rhs1 + tau*(K*rhoq));
        rhs2 = -invMVrT*(rhs2 + tau*(K*rhouq));
        rhs3 = -invMVrT*(rhs3 + tau*(K*rhovq));
        rhs4 = -invMVrT*(rhs4 + tau*(K*Eq));
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        res4 = rk4a(INTRK)*res4 + dt*rhs4;
        
        rho  = rho  + rk4b(INTRK)*res1;
        rhou = rhou + rk4b(INTRK)*res2;
        rhov = rhov + rk4b(INTRK)*res3;
        E    = E    + rk4b(INTRK)*res4;        
    end
        
    
    if (mod(i,5)==0 || i==Nsteps)
        surf(x,y,reshape(Vrp*rhou,Np1D,Np1D))
        view(2)
        colorbar
        shading interp
        title(sprintf('step = %d / %d\n',i,Nsteps))
        drawnow
    end
end

err = Vrp*rhou-Usnap((1:Np)+Np,end);
% surf(x,y,reshape(err,2*N+1,2*N+1));shading interp
fprintf('final err = %g\n',sqrt(sum(dx^2*(err).^2))/sqrt(sum(dx^2*(Vrp*rhou).^2)))

