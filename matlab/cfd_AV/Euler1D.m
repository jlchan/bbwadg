clear
clear -global
Globals1D;

CFL = .125;
N = 2;
K1D = 64;
plotErrEst = 0;

FinalTime = .5;
opt = 2;

global tau
tau = 1;

r = JacobiGL(0,0,N);
% r = JacobiGQ(0,0,N);

[rq wq] = JacobiGQ(0,0,N); rq = rq*(1-1e-11); 


% [rq wq] = JacobiGQ(0,0,N+1); rq = rq*(1-1e-11); 
% [rq wq] = JacobiGL(0,0,N+4); rq = rq*(1-1e-11);
% [rq wq] = JacobiGQ(0,0,N+1);


if opt==1 || opt==0
    [Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
elseif opt==2
    [Nv, VX, K, EToV] = MeshGen1D(-.5,.5,K1D);
    FinalTime = .2;
end

% Initialize solver and construct grid and metric
StartUp1D;

global Ef Vq Pq Lf Vf wq E
V = Vandermonde1D(N,r);
Dr = GradVandermonde1D(N,r)/V;
Vq = Vandermonde1D(N,rq)/V;
Vf = Vandermonde1D(N,[-1 1])/V;
xq = Vq*x;

rp = linspace(-1,1,50); 
Vp = Vandermonde1D(N,rp)/V;
xp = Vandermonde1D(N,rp)/Vandermonde1D(N,JacobiGL(0,0,N))*x;

M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq));
VqPq = Vq*Pq;
wf = [1;1];
Lf = Vq*(M\(Vf'*diag(wf)));

global nxJ rxJ
rxJ = rx.*J;
nxJ = nx;

W = diag([wq;wf]);

B = diag([-1,1]);
E = Vf*Pq;
Qr = Pq'*M*Dr*Pq;
global QNr
QNr = .5*[Qr-Qr' E'*B;
    -B*E B];

QNr = .5*(QNr-QNr'); %skew form

global VPN
VNP = [Vq;Vf]*Pq;
VPN = Vq*(M\[Vq;Vf]');

%% extra dissipation (subcell?)

[rq2 wq2] = JacobiGQ(0,0,N+2); 
Vq2 = Vandermonde1D(N,rq2)/V;
Pq2 = (Vq2'*diag(wq2)*Vq2)\(Vq2'*diag(wq2));
xq2 = Vq2*x;

%% maps

global mapP mapB
mapM = reshape(1:Nfp*Nfaces*K,Nfp*Nfaces,K);
mapP = mapM;
for e = 1:K
    mapP(:,e) = [(Nfp*Nfaces)*e-2, Nfp*Nfaces*(e+1)-1];
end
mapP(1) = 1; mapP(end) = Nfp*Nfaces*K;
mapP = reshape(mapP,Nfp*Nfaces,K);
mapB = [1; Nfp*Nfaces*K];


%% fluxes

global gamma
gamma = 1.4;
rhoe = @(rho,m,E) E - .5*m.^2./rho;
s = @(rho,m,E) log((gamma-1).*rhoe(rho,m,E)./(rho.^gamma));

global V1 V2 V3
V1 = @(rho,m,E) (-E + rhoe(rho,m,E).*(gamma + 1 - s(rho,m,E)))./(rhoe(rho,m,E));
V2 = @(rho,m,E) (m)./(rhoe(rho,m,E));
V3 = @(rho,m,E) (-rho)./(rhoe(rho,m,E));

sV = @(V1,V2,V3) gamma - V1 + V2.^2./(2*V3);
rhoeV  = @(V1,V2,V3) ((gamma-1)./((-V3).^gamma)).^(1/(gamma-1)).*exp(-sV(V1,V2,V3)/(gamma-1));
U1 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(-V3);
U2 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(V2);
U3 = @(V1,V2,V3) rhoeV(V1,V2,V3).*(1-V2.^2./(2*V3));

global avg pfun beta
avg = @(uL,uR) .5*(uL+uR);
pfun = @(rho,u,E) (gamma-1)*(E-.5*rho.*u.^2);
beta = @(rho,u,E) rho./(2*pfun(rho,u,E));

global f1 f2 f3
pavg = @(rhoL,rhoR,uL,uR,EL,ER) avg(rhoL,rhoR)./(2*avg(beta(rhoL,uL,EL),beta(rhoR,uR,ER)));
f1 = @(rhoL,rhoR,uL,uR,EL,ER) logmean(rhoL,rhoR).*avg(uL,uR);
f2 = @(rhoL,rhoR,uL,uR,EL,ER) pavg(rhoL,rhoR,uL,uR,EL,ER) + avg(uL,uR).*f1(rhoL,rhoR,uL,uR,EL,ER);
f3 = @(rhoL,rhoR,uL,uR,EL,ER) f1(rhoL,rhoR,uL,uR,EL,ER)...
    .*(1./(2*(gamma-1).*logmean(beta(rhoL,uL,EL),beta(rhoR,uR,ER))) - .5*avg(uL.^2,uR.^2)) ...
    + avg(uL,uR).*f2(rhoL,rhoR,uL,uR,EL,ER);

% % ranocha EC correction
% f1 = @(rhoL,rhoR,uL,uR,EL,ER) logmean(rhoL,rhoR).*avg(uL,uR);
% f2 = @(rhoL,rhoR,uL,uR,EL,ER) avg(pfun(rhoL,uL,EL),pfun(rhoR,uR,ER)) + avg(uL,uR).*f1(rhoL,rhoR,uL,uR,EL,ER);
% f3 = @(rhoL,rhoR,uL,uR,EL,ER) f1(rhoL,rhoR,uL,uR,EL,ER)...
%     .*(1./(2*(gamma-1).*logmean(beta(rhoL,uL,EL),beta(rhoR,uR,ER))) - .5*avg(uL.^2,uR.^2)) ...
%     + avg(uL,uR).*f2(rhoL,rhoR,uL,uR,EL,ER) + avg(pfun(rhoL,uL,EL),pfun(rhoR,uR,ER)).*avg(uL,uR) - (pfun(rhoL,uL,EL)-pfun(rhoR,uR,ER)).*(uL-uR)/4;

%%
res1 = 0;
res2 = 0;
res3 = 0;

CN = (N+1)^2/2;
h = 2/K1D;
dt = CFL * h/CN;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

global rhoL rhoR mL mR EL ER
global useBC; 

if opt==0
    
    % sine solution
    t = 0;
    rhoex = @(x) (2 + sin(pi*(x - t)));
    uex = @(x) ones(size(x));
    pex = @(x) ones(size(x));
    mex = @(x) rhoex(x).*uex(x);
    
    Eex = @(x) pex(x)./(gamma-1) + .5*rhoex(x).*uex(x).^2;
    
    rhoq = rhoex(xq);
    uq = uex(xq);
    pq = pex(xq);
    mq = mex(xq);
    Eq = Eex(xq);
    
    useBC = 0;
    mapP(1) = mapM(end); mapP(end) = mapM(1);
    vmapP(1) = vmapM(end); vmapP(end) = vmapM(1);
    
elseif opt==1
    
    % pulse condition
    x0 = 0;
    dU = -.5*sin(pi*xq); 
    dU = (abs(xq-x0) < .5);
    rhoq = (2 + dU);
    pq = rhoq.^gamma;
    uq = 0*dU;
    Eq = pq/(gamma-1) + .5*rhoq.*uq.^2;
    mq = rhoq.*uq;
    
    useBC = 0;
    mapP(1) = mapM(end); mapP(end) = mapM(1);
    vmapP(1) = vmapM(end); vmapP(end) = vmapM(1);
    
elseif opt==2
    
    % BCs
    useBC = 1;
    rhoL = 1; rhoR = .125;
    pL = 1; pR = .1;
    %     rhoR = .1; pR = .05;
    uL = 0; uR = 0;
    mL = 0; mR = 0;
    EL = pL/(gamma-1) + .5*rhoL*uL.^2;
    ER = pR/(gamma-1) + .5*rhoR*uR.^2;
    
    % sod initial condition
    rhoq = (rhoL*(xq < 0) + rhoR*(xq > 0));
    pq = (pL*(xq < 0) + pR*(xq > 0));
    uq = 0*xq;
    Eq = pq/(gamma-1) + .5*rhoq.*uq.^2;
    mq = rhoq.*uq;
    
end

rho = VqPq*(rhoq);
u = VqPq*(uq);
p = VqPq*(pq);
E = VqPq*(pq/(gamma-1) + .5*rhoq.*uq.^2);
m = VqPq*(rhoq.*uq);

wJq = diag(wq)*(Vq*J);

figure
for i = 1:Nsteps
    for INTRK = 1:5                
        
        % project entropy variables
        q1 = VNP*V1(rho,m,E);
        q2 = VNP*V2(rho,m,E);
        q3 = VNP*V3(rho,m,E);
        
        rhoN = U1(q1,q2,q3);
        mN   = U2(q1,q2,q3);
        EN   = U3(q1,q2,q3);
        uN   = mN./rhoN;
        
        [rhs1 rhs2 rhs3] = rhsEuler(rhoN,uN,EN);
        
        if INTRK==5
            v1 = VqPq*V1(rho,m,E);
            v2 = VqPq*V2(rho,m,E);
            v3 = VqPq*V3(rho,m,E);
            rhstest = sum(sum(wJq.*(v1.*rhs1 + v2.*rhs2 + v3.*rhs3)));            
        end                
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        rho = rho + rk4b(INTRK)*res1;
        m   = m   + rk4b(INTRK)*res2;
        E   = E   + rk4b(INTRK)*res3;                
    end
        
    Sq = -rho.*s(rho,m,E);
    S(i) = sum(sum(wJq.*Sq));    
    
    if mod(i,5)==0 || i==Nsteps
        
        u = m./rho;
        p = (gamma-1)*(E-.5*rho.*u.^2);
        
        rhohat = V\(Pq*rho);
        mhat   = V\(Pq*m);
        Ehat   = V\(Pq*E);
        modalEst = (rhohat(end,:).^2 + mhat(end,:).^2 + Ehat(end,:).^2)./sum(rhohat.^2+mhat.^2+Ehat.^2,1);        
        modalEst = repmat(modalEst,length(rp),1);                
        
        % dissipation difference terms: coarse/fine
        v1c = Vq2*Pq*V1(rho,m,E);
        v2c = Vq2*Pq*V2(rho,m,E);
        v3c = Vq2*Pq*V3(rho,m,E);
        v1f = Vq2*Pq2*V1(Vq2*Pq*rho,Vq2*Pq*m,Vq2*Pq*E);
        v2f = Vq2*Pq2*V2(Vq2*Pq*rho,Vq2*Pq*m,Vq2*Pq*E);
        v3f = Vq2*Pq2*V3(Vq2*Pq*rho,Vq2*Pq*m,Vq2*Pq*E);
        dv1 = v1f-v1c;
        dv2 = v2f-v2c;
        dv3 = v3f-v3c;
        
        % eval on fine quadrature
        u1c = U1(v1c,v2c,v3c);
        u2c = U2(v1c,v2c,v3c);
        u3c = U3(v1c,v2c,v3c);
        u1f = U1(v1f,v2f,v3f);
        u2f = U2(v1f,v2f,v3f);
        u3f = U3(v1f,v2f,v3f);
        du1 = u1f-u1c;
        du2 = u2f-u2c;
        du3 = u3f-u3c;
        
        wJq2 = diag(wq2)*(Vq2*J);        
        dVdU = (dv1.*du1 + dv2.*du2 + dv3.*du3);        
        
        % How to show dS/dU*U > 0?  
        Savg = abs(sum(wJq.*(V1(rho,m,E).*rho + V2(rho,m,E).*m + V3(rho,m,E).*E),1));
        dVdUEst = repmat(abs(sum(wJq2.*dVdU,1))./Savg,length(rp),1);                
        
        v1f = Vf*Pq*V1(rho,m,E);
        v2f = Vf*Pq*V2(rho,m,E);
        v3f = Vf*Pq*V3(rho,m,E);
        dv1f = V1(Vf*Pq*rho,Vf*Pq*m,Vf*Pq*E)-v1f;
        dv2f = V2(Vf*Pq*rho,Vf*Pq*m,Vf*Pq*E)-v2f;
        dv3f = V3(Vf*Pq*rho,Vf*Pq*m,Vf*Pq*E)-v3f;
        du1f = Vf*Pq*rho - U1(v1f,v2f,v3f);
        du2f = Vf*Pq*m - U2(v1f,v2f,v3f);
        du3f = Vf*Pq*E - U3(v1f,v2f,v3f);
        %         du1_Lift = Lf*(Vf*Pq*rho - U1(v1f,v2f,v3f));
        %         du2_Lift = Lf*(Vf*Pq*m - U2(v1f,v2f,v3f));
        %         du3_Lift = Lf*(Vf*Pq*E - U3(v1f,v2f,v3f));
        %         dVdULift = (V1(rho,m,E)).*du1_Lift + (V2(rho,m,E)).*du2_Lift + (V3(rho,m,E)).*du3_Lift;
        
        dVdULift = Lf*(dv1f.*du1f + dv2f.*du2f + dv3f.*du3f);
%         dVdULift = (Vq*Pq2*dv1).*(Lf*du1f) + (Vq*Pq2*dv2).*(Lf*du2f) + (Vq*Pq2*dv3).*(Lf*du3f);
        dVdULiftEst = repmat(abs(sum(wJq.*dVdULift,1))./Savg,length(rp),1);
                
        dULift = CN*sum(du1f.^2 + du2f.^2 + du3f.^2);
        dULiftEst = repmat(abs(sum(dULift,1))./sum(wJq.*(rho.^2 + m.^2 + E.^2)),length(rp),1);
        
        if plotErrEst
            t = i*dt;
            tol = -10;           
            
            subplot(2,1,1)            
            vv = modalEst;             
            vv = log10(vv);
            vv(vv<tol) = nan;
            color_line3(xp,ones(size(xp))*t,vv,vv,'.')
            title('Modal estimator')
            
%             subplot(3,1,2)
%             vv = abs(dVdUEst); 
%             vv = vv/max(vvref(:)); 
%             vv(abs(vv)<tol) = nan;
%             color_line3(xp,ones(size(xp))*t,vv,vv,'.')
%             title('Volume entropy est')
            
            subplot(2,1,2)            
            vv = abs(dULiftEst);
            vv = log10(vv);
            vv(vv<tol) = nan;
            color_line3(xp,ones(size(xp))*t,vv,vv,'.')
            title('Lifted entropy est')
        else
            
%             rhop = Vp*Pq*rho;
%             mp = Vp*Pq*m;
%             Ep = Vp*Pq*E;
Savg = repmat(.5*wq'*s(rho,m,E),length(rp),1);
            plot(xp,Savg,'b-','linewidth',2)
            hold on
            plot(xp,Vp*Pq*rho,'k-','linewidth',2)
            
%             plot(xq,VqPq*u,'r-','linewidth',2)
%             plot(xq,VqPq*p,'k-','linewidth',2)
                        
            plot(xp,modalEst/max(abs(modalEst(:))),'ro-')
%             vv = abs(dVdULiftEst);
%             plot(xp,vv/max(vv(:)),'b.-')
%             %         plot(xq,(Vq*Pq2*dVAdU)/magAV,'r-')            
            title(sprintf('Time = %f, rhstest = %g',dt*i,rhstest))
        
        end
        
                
        %         axis([-5 5 -1 7])
        hold off
        drawnow limitrate
    end
end

if opt==0
    % sine solution
    t = FinalTime;
    rhoex = @(x) (2 + sin(pi*(x - t)));
    uex = @(x) ones(size(x));
    pex = @(x) ones(size(x));
    mex = @(x) rhoex(x).*uex(x);
    Eex = @(x) pex(x)./(gamma-1) + .5*rhoex(x).*uex(x).^2;
    
    e1 = rhoex(xq) - Vq*rho;
    e2 = mex(xq) - Vq*m;
    e3 = Eex(xq) - Vq*E;
    L2err = sqrt(sum(sum(wJq.*(e1.^2+e2.^2+e3.^2))));
    L2err
end


function [rhs1 rhs2 rhs3] = rhsEuler(rhoN,uN,EN)

Globals1D
global rhoL rhoR mL mR EL ER gamma
global mapP mapB
global f1 f2 f3
global QNr Pq Lq
global Vf wq Vq

rhoM = rhoN(end-1:end,:);
uM = uN(end-1:end,:);
EM = EN(end-1:end,:);
mM = rhoM.*uM;
rhoP = reshape(rhoM(mapP),Nfp*Nfaces,K);
uP = reshape(uM(mapP),Nfp*Nfaces,K);
EP = reshape(EM(mapP),Nfp*Nfaces,K);
mP = rhoP.*uP;
pM = (gamma-1)*(EM-.5*rhoM.*uM.^2);

global useBC
if useBC
    rhoP(mapB) = [rhoL rhoR];
    uP(mapB) = [mL/rhoL mR/rhoR];
    mP(mapB) = [mL mR];
    EP(mapB) = [EL ER];
end

global tau

% local lax penalty
lm   = (abs(uM) + sqrt(gamma*pM./rhoM));
Lfc  = max(lm(mapP),lm);

% compute fluxes
global nxJ PN Lf
f1f = nxJ.*f1(rhoM,rhoP,uM,uP,EM,EP) - .5*tau*Lfc.*((rhoP-rhoM));
f2f = nxJ.*f2(rhoM,rhoP,uM,uP,EM,EP) - .5*tau*Lfc.*((mP-mM));
f3f = nxJ.*f3(rhoM,rhoP,uM,uP,EM,EP) - .5*tau*Lfc.*((EP-EM));

r1 = zeros(size(QNr,1),K);
r2 = zeros(size(QNr,1),K);
r3 = zeros(size(QNr,1),K);
for e = 1:K
    [rhox rhoy] = meshgrid(rhoN(:,e));
    [ux uy] = meshgrid(uN(:,e));
    [Ex Ey] = meshgrid(EN(:,e));
    
    FS1 = f1(rhox,rhoy,ux,uy,Ex,Ey);
    FS2 = f2(rhox,rhoy,ux,uy,Ex,Ey);
    FS3 = f3(rhox,rhoy,ux,uy,Ex,Ey);
    r1(:,e) = 2*sum(QNr.*FS1,2); % rxJ = 1 in 1D
    r2(:,e) = 2*sum(QNr.*FS2,2);
    r3(:,e) = 2*sum(QNr.*FS3,2);
end

global VPN
rhs1 = VPN*r1 + Lf*f1f;
rhs2 = VPN*r2 + Lf*f2f;
rhs3 = VPN*r3 + Lf*f3f;

rhs1 = -rhs1./J(1);
rhs2 = -rhs2./J(1);
rhs3 = -rhs3./J(1);


end