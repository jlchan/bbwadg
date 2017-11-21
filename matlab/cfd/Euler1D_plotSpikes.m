clear
Globals1D;

projectV = 1;
CFL = .125/4;
N = 2;
K1D = 40;
FinalTime = 1;
% FinalTime = 3.e-04;
opt = 3;

global tau
tau = 1;

r = JacobiGL(0,0,N);
r = JacobiGQ(0,0,N);

[rq wq] = JacobiGL(0,0,N); rq = .999999999*rq;
[rq wq] = JacobiGQ(0,0,N+1);


% % include boundary nodes for extraction
rq = [-1;rq;1];
wq = [0;wq;0];
Nq = length(rq);

% evaluate error
rq2 = rq; wq2 = wq;
% [rq2 wq2] = JacobiGQ(0,0,N+4);
% [rq2 wq2] = JacobiGL(0,0,N);

if opt==1
    [Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
elseif opt==2
    [Nv, VX, K, EToV] = MeshGen1D(-.5,.5,K1D);
elseif opt==3
    [Nv, VX, K, EToV] = MeshGen1D(-5,5,K1D);
end

% Initialize solver and construct grid and metric
StartUp1D;

global mapP mapB
mapM = reshape(1:Nfp*Nfaces*K,Nfp*Nfaces,K);
mapP = mapM;
for e = 1:K
    mapP(:,e) = [(Nfp*Nfaces)*e-2, Nfp*Nfaces*(e+1)-1];
end
mapP(1) = 1; mapP(end) = Nfp*Nfaces*K;
mapP = reshape(mapP,Nfp*Nfaces,K);
mapB = [1; Nfp*Nfaces*K];

global VqPq tf Dq Vq Pq Dfq
V = Vandermonde1D(N,r);
Dr = GradVandermonde1D(N,r)/V;
Vq = Vandermonde1D(N,rq)/V;
Vf = Vandermonde1D(N,[-1 1])/V;
Vq2 = Vandermonde1D(N,rq2)/V;

M = Vq'*diag(wq)*Vq;
LIFT = M\Vf';

Pq = M\(Vq'*diag(wq));
xq =  Vandermonde1D(N,rq)/Vandermonde1D(N,JacobiGL(0,0,N))*x;
Pq(abs(Pq)<1e-8) = 0;

VqPq = Vq*Pq;
tf = zeros(2,length(rq)); tf(1) = 1; tf(end) = 1;
Dq = Vq*Dr*Pq;
Dfq = .5*Vq*LIFT*diag([-1,1])*(tf - tf*VqPq);

rp = linspace(-1,1,50); F = eye(N+1);
Vp = Vandermonde1D(N,rp)/V;
xp = Vandermonde1D(N,rp)/Vandermonde1D(N,JacobiGL(0,0,N))*x;

rp0 = 0; F0 = V*diag([1;zeros(N,1)])/V;
Vp0 = Vandermonde1D(N,rp0)/V;
xp0 = Vandermonde1D(N,rp0)/Vandermonde1D(N,JacobiGL(0,0,N))*x;


%% fluxes

global gamma
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

global f1 f2 f3
f1 = @(rhoL,rhoR,uL,uR,EL,ER) logmean(rhoL,rhoR).*avg(uL,uR);
f2 = @(rhoL,rhoR,uL,uR,EL,ER) avg(rhoL,rhoR)./(2*avg(beta(rhoL,uL,EL),beta(rhoR,uR,ER))) + avg(uL,uR).*f1(rhoL,rhoR,uL,uR,EL,ER);
f3 = @(rhoL,rhoR,uL,uR,EL,ER) f1(rhoL,rhoR,uL,uR,EL,ER)...
    .*(1./(2*(gamma-1).*logmean(beta(rhoL,uL,EL),beta(rhoR,uR,ER))) - .5*avg(uL.^2,uR.^2)) ...
    + avg(uL,uR).*f2(rhoL,rhoR,uL,uR,EL,ER);


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
global useBC; useBC = 1;

if opt==1
    % pulse condition
    x0 = 0;
    rho = Pq*(2 + (abs(xq-x0) < .5));
    p = rho.^gamma;
    u = 0*rho;
    E = p/(gamma-1) + .5*rho.*u.^2;
    m = rho.*u;
    
    useBC = 0;
    mapP(1) = mapM(end); mapP(end) = mapM(1);
    vmapP(1) = vmapM(end); vmapP(end) = vmapM(1);
    
elseif opt==2
    % sod initial condition
    rho = Pq*(1*(xq < 0) + .125*(xq > 0));
    p = Pq*(1*(xq < 0) + .1*(xq > 0));
    u = 0*rho;
    E = p/(gamma-1) + .5*rho.*u.^2;
    m = rho.*u;
    
    % BCs
    rhoL = 1; rhoR = .125;
    pL = 1; pR = .1;
    uL = 0; uR = 0;
    mL = 0; mR = 0;
    EL = pL/(gamma-1) + .5*rhoL*uL.^2;
    ER = pR/(gamma-1) + .5*rhoR*uR.^2;
    
elseif opt==3
    
    rhoRex = @(x) 1 + .2*sin(5*x);
    rhoL = 3.857143; rhoR = rhoRex(x(end));
    uL   = 2.629369; uR = 0;
    pL   = 10.3333;  pR = 1;
    
%     rhoRex = @(x) 1 + .2*sin(5*x);
%     rhoL = 3.857143; rhoR = rhoRex(x(end));
%     uL   = 2.629369; uR = 0;
%     pL   = 2;  pR = 1;
        
%     useBC = 0;
%     mapP(1) = mapM(end); mapP(end) = mapM(1);
%     vmapP(1) = vmapM(end); vmapP(end) = vmapM(1);
    
    rhoq = rhoL*(xq < -4) + (rhoRex(xq)).*(xq >= -4);
    uq = uL*(xq < -4);
    pq = pL*(xq < -4) + pR*(xq >= -4);
    rho = Pq*(rhoq);
    u = Pq*(uq);
    p = Pq*(pq);
    E = Pq*(pq/(gamma-1) + .5*rhoq.*uq.^2);
    m = Pq*(rhoq.*uq);
    
    mL = rhoL*uL; mR = rhoR*uR;
    EL = pL/(gamma-1) + .5*rhoL*uL.^2;
    ER = pR/(gamma-1) + .5*rhoR*uR.^2;
end

wJq = diag(wq)*(Vq*J);

figure(1)
for i = 1:Nsteps    
    
    for INTRK = 1:5
        
        % interpolate to quadrature
        rhoq = Vq*rho;
        mq = Vq*m;
        Eq = Vq*E;
        
        % project entropy variables
        q1 = VqPq*V1(rhoq,mq,Eq);
        q2 = VqPq*V2(rhoq,mq,Eq);
        q3 = VqPq*V3(rhoq,mq,Eq);
        
        if projectV
            rhoqtmp = U1(q1,q2,q3);
            mqtmp   = U2(q1,q2,q3);
            Eqtmp   = U3(q1,q2,q3);

%             % apply limiting to projected variables: 

% rhoavg = sum(wq.*rhoqtmp,1)/sum(wq);            
% t = @(u,p,pavg) repmat(min(1,min(abs((max(u,[],1)-pavg)./(max(p,[],1) - pavg)), abs((min(u,[],1)-pavg)./(min(u,[],1) - pavg)))),Nq,1);
% rhoq = repmat(rhoavg,Nq,1) + t(rhoq,rhoqtmp,rhoavg).*(rhoqtmp - repmat(rhoavg,Nq,1));

            
%             t = min(abs((max(rhoq,[],1)-rhoavg)./(max(rhoqtmp,[],1) - rhoavg)), abs((min(rhoq,[],1)-rhoavg)./(min(rhoqtmp,[],1) - rhoavg)));
%             t = repmat(min(t,1),Nq,1);
%             rhoavg = repmat(rhoavg,Nq,1);            
%             rhoq = rhoavg + t.*(rhoqtmp-rhoavg);
%             
%             mavg = sum(wq.*mqtmp,1)/sum(wq);            
%             t = min(abs((max(mq,[],1)-mavg)./(max(mqtmp,[],1) - mavg)), abs((min(mq,[],1)-mavg)./(min(mqtmp,[],1) - mavg)));
%             t = repmat(min(t,1),Nq,1);
%             mavg = repmat(mavg,Nq,1);
%             mq = mavg + t.*(mqtmp-mavg);
%             
%             Eavg = sum(wq.*Eqtmp,1)/sum(wq);
%             t = min(abs((max(Eq,[],1)-Eavg)./(max(Eqtmp,[],1) - Eavg)), abs((min(Eq,[],1)-Eavg)./(min(Eqtmp,[],1) - Eavg)));
%             t = repmat(min(t,1),Nq,1);
%             Eavg = repmat(Eavg,Nq,1);
%             Eq = Eavg + t.*(Eqtmp-Eavg);
            
        end                
        
        [rhs1 rhs2 rhs3] = rhsEuler(rhoq,mq,Eq,INTRK);
        
        if (INTRK==5)
            rhstest(i) = sum(sum(wJq.*(q1.*(Vq*rhs1) + q2.*(Vq*rhs2) + q3.*(Vq*rhs3))));
        end
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        rho = rho + rk4b(INTRK)*res1;
        m = m + rk4b(INTRK)*res2;
        E = E + rk4b(INTRK)*res3;
        
        if 1
            rhop = Vp*F*rho;
            mp = Vp*F*m;
            Ep = Vp*F*E;
            up = mp./rhop;
            pp = (gamma-1)*(Ep-.5*rhop.*up.^2);
            plot(xp,rhop,'-','linewidth',2)
            hold on
            plot(xp,up,'--','linewidth',2)
            plot(xp,pp,'--','linewidth',2)
            title(sprintf('Time = %f, tstep = %d, stage = %d',dt*i,i,INTRK))
            %         axis([-5 5 -1 7])            
            hold off
            drawnow 
            if INTRK==5 && i==4
                
                % project entropy variables                
                rhoq = Vq*rho;
                mq = Vq*m;
                Eq = Vq*E;                                  
                v1 = V1(rhoq,mq,Eq);
                v2 = V2(rhoq,mq,Eq);
                v3 = V3(rhoq,mq,Eq);
                q1 = VqPq*v1;
                q2 = VqPq*v2;
                q3 = VqPq*v3;
                
                q3avg = repmat(wq'*q3,Nq,1)/sum(wq);
                q3 = q3avg + .8*(q3 - q3avg);
                
%                 plot(xq,q3,'o--')
%                 hold on
%                 plot(xq,q3lim,'x--')
%                 keyboard
                
                rhoqproj = U1(q1,q2,q3);
                mqproj   = U2(q1,q2,q3);
                Eqproj   = U3(q1,q2,q3);
                
                % plotting entropy variables                
                rhop = Vp*rho;
                mp = Vp*m;
                Ep = Vp*E;
                v1p = V1(rhop,mp,Ep);
                v2p = V2(rhop,mp,Ep);
                v3p = V3(rhop,mp,Ep);
                rhoqprojp = U1(Vp*Pq*q1,Vp*Pq*q2,Vp*Pq*q3);
                mqprojp   = U2(Vp*Pq*q1,Vp*Pq*q2,Vp*Pq*q3);
                Eqprojp   = U3(Vp*Pq*q1,Vp*Pq*q2,Vp*Pq*q3);
                
                up = mp./rhop;
                pp = (gamma-1)*(Ep - .5*rhop.*up.^2);
                uprojp = mqprojp./rhoqprojp;
                pprojp = (gamma-1)*(Eqprojp - .5*rhoqprojp.*uprojp.^2);
                
                h1 = plot(xq,rhoq,'bo-','linewidth',2);
                hold on
                h2 = plot(xq,rhoqproj,'rx-','linewidth',2);
                
                keyboard
                clf
                h1 = plot(xp,v1p,'b-','linewidth',2);
                hold on
                h2 = plot(xp,Vp*Pq*q1,'r--','linewidth',2);
                legend([h1(1) h2(1)],{'v_1','Projected v_1'},'FontSize',15,'Location','NorthEast')
                set(gca,'fontsize',15);grid on
                axis([-4.75 -1.75 0 2.5])                
%                 print(gcf,'-dpng','~/Desktop/bbwadg/docs/splitform/figs/sineShockQ1Compare.png')
                
                clf
                h1 = plot(xp,v2p,'b-','linewidth',2);
                hold on
                h2 = plot(xp,Vp*Pq*q2,'r--','linewidth',2);
                legend([h1(1) h2(1)],{'v_2','Projected v_2'},'FontSize',15,'Location','Best')
                set(gca,'fontsize',15);grid on
                axis([-4.75 -1.75 -.25 .5])                
%                 print(gcf,'-dpng','~/Desktop/bbwadg/docs/splitform/figs/sineShockQ2Compare.png')
                
                clf
                h1 = plot(xp,v3p,'b-','linewidth',2);
                hold on
                h2 = plot(xp,Vp*Pq*q3,'r--','linewidth',2);
                legend([h1(1) h2(1)],{'v_3','Projected v_3'},'FontSize',15,'Location','Best')
                set(gca,'fontsize',15);grid on
                axis([-4.75 -1.75 -.65 0])                
%                 print(gcf,'-dpng','~/Desktop/bbwadg/docs/splitform/figs/sineShockQ3Compare.png')
                
                clf
                h1 = plot(xp,rhop,'b-','linewidth',2);
                hold on
                h2 = plot(xp,rhoqprojp,'r--','linewidth',2);
                axis([-4.75 -1.75 0 22])                
                legend([h1(1) h2(1)],{'Density (original)','Density (projected entropy variables)'},'FontSize',15,'Location','NorthEast')
                set(gca,'fontsize',15);grid on
%                 print(gcf,'-dpng','~/Desktop/bbwadg/docs/splitform/figs/sineShockDensityCompare.png')
                
                clf
                h1 = plot(xp,up,'b-','linewidth',2);
                hold on
                h2 = plot(xp,uprojp,'r--','linewidth',2);
                axis([-4.75 -1.75 -1 5])                
                legend([h1(1) h2(1)],{'Velocity (original)','Velocity (projected entropy variables)'},'FontSize',15,'Location','NorthEast')
                set(gca,'fontsize',15);grid on
                print(gcf,'-dpng','~/Desktop/bbwadg/docs/splitform/figs/sineShockVelocityCompare.png')
                
                clf
                h1 = plot(xp,pp,'b-','linewidth',2);
                hold on
                h2 = plot(xp,pprojp,'r--','linewidth',2);
                axis([-4.75 -1.75 -10 100])                
                legend([h1(1) h2(1)],{'Pressure (original)','Pressure (projected entropy variables)'},'FontSize',15,'Location','NorthEast')
                set(gca,'fontsize',15);grid on
                print(gcf,'-dpng','~/Desktop/bbwadg/docs/splitform/figs/sineShockPressureCompare.png')
                
                keyboard

            end
        end
    end
    
    
    
    if mod(i,10)==0 || i==Nsteps
            rhop = Vp*F*rho;
            mp = Vp*F*m;
            Ep = Vp*F*E;
            up = mp./rhop;
            pp = (gamma-1)*(Ep-.5*rhop.*up.^2);
            plot(xp,rhop,'-','linewidth',2)
            hold on
%             plot(xp,up,'--','linewidth',2)
            plot(xp,pp,'--','linewidth',2)
            title(sprintf('Time = %f',dt*i))
            %         axis([-5 5 -1 7])            
            hold off
            drawnow            
        end
    
end

return
%%
rhop = Vp*F*rho;
mp = Vp*F*m;
Ep = Vp*F*E;
up = mp./rhop;
pp = (gamma-1)*(Ep-.5*rhop.*up.^2);

xp = [xp;nan(1,size(rhop,2))];
rhop = [rhop;nan(1,size(rhop,2))];
up = [up;nan(1,size(rhop,2))];
pp = [pp;nan(1,size(rhop,2))];


clf
plot(xp(:),rhop(:),'b-','linewidth',3)
hold on
plot(xp(:),pp(:),'r--','linewidth',3)
plot(xp0,Vp0*F0*rho,'ko','linewidth',3,'markersize',14,'MarkerFaceColor',[.49 1 .63])

rho0 = Vp0*F0*rho;
E0 = Vp0*F0*E;
m0 = Vp0*F0*m;
p0 = (gamma-1)*(E0 - .5*m0.^2./rho0);
plot(xp0,p0,'ko','linewidth',3,'markersize',14,'MarkerFaceColor',[.49 1 .63])
% plot(xp0,Vp0*F0*rho,'kx','linewidth',2)
hold on
% plot(xp(:),pp(:),'r--','linewidth',2)
% axis([-4.2 -3.8 -1 12])
axis([-.5 .5 -.1 1.1])
legend({'Density','Pressure'},'FontSize',32)
grid on; set(gca,'fontsize',32)
% print(gcf,'-dpng','~/Desktop/bbwadg/docs/splitform/figs/sodGLL.png')
%print(gcf,'-dpng','~/Desktop/bbwadg/docs/splitform/figs/sineshockGLL.png')
% print(gcf,'-dpng','~/Desktop/bbwadg/docs/splitform/figs/sineshockGQN2.png')

%%
return



figure(2)
semilogy(dt*(1:Nsteps),rhstest,'x')

function [rhs1 rhs2 rhs3] = rhsEuler(rhoq,mq,Eq,step)

Globals1D
global rhoL rhoR mL mR EL ER gamma
global mapP mapB
global VqPq tf Dq Vq Pq Dfq
global f1 f2 f3

uq = mq./rhoq;
pq = (gamma-1)*(Eq-.5*rhoq.*uq.^2);

Nq = size(rhoq,1);
rhoM = rhoq([1 Nq],:); uM = uq([1 Nq],:); EM = Eq([1 Nq],:);
mM = rhoM.*uM;
rhoP = reshape(rhoM(mapP),Nfp*Nfaces,K);
uP = reshape(uM(mapP),Nfp*Nfaces,K);
EP = reshape(EM(mapP),Nfp*Nfaces,K);
mP = rhoP.*uP;
pM = pq([1 Nq],:);

global useBC
if useBC
    rhoP(mapB) = [rhoL rhoR];
    uP(mapB) = [mL/rhoL mR/rhoR];
    mP(mapB) = [mL mR];
    EP(mapB) = [EL ER];
end

% compute fluxes
f1f = f1(rhoM,rhoP,uM,uP,EM,EP);
f2f = f2(rhoM,rhoP,uM,uP,EM,EP);
f3f = f3(rhoM,rhoP,uM,uP,EM,EP);

rhs1 = zeros(Np,K);
rhs2 = zeros(Np,K);
rhs3 = zeros(Np,K);
for e = 1:K
    [rhox rhoy] = meshgrid(rhoq(:,e));
    [ux uy] = meshgrid(uq(:,e));
    [Ex Ey] = meshgrid(Eq(:,e));
    
    [rhoxf rhoyf] = meshgrid(rhoq(:,e),tf*rhoq(:,e));
    [uxf uyf] = meshgrid(uq(:,e),tf*uq(:,e));
    [Exf Eyf] = meshgrid(Eq(:,e),tf*Eq(:,e));
    
    flux1 = .5*diag([-1 1])*f1f(:,e);
    flux2 = sum((.5*diag([-1 1])*(-tf*VqPq)).*f1(rhoxf,rhoyf,uxf,uyf,Exf,Eyf),2);
    rhs1(:,e) = rx(:,e).*(Pq*sum(Dq.*f1(rhox,rhoy,ux,uy,Ex,Ey),2)) + Fscale(1,e)*(Pq*sum(Dfq.*f1(rhox,rhoy,ux,uy,Ex,Ey),2)) + LIFT*(Fscale(:,e).*(flux1+flux2));
    
    flux1 = .5*diag([-1 1])*f2f(:,e);
    flux2 = sum((.5*diag([-1 1])*(-tf*VqPq)).*f2(rhoxf,rhoyf,uxf,uyf,Exf,Eyf),2);
    rhs2(:,e) = rx(:,e).*(Pq*sum(Dq.*f2(rhox,rhoy,ux,uy,Ex,Ey),2)) + Fscale(1,e)*(Pq*sum(Dfq.*f2(rhox,rhoy,ux,uy,Ex,Ey),2)) + LIFT*(Fscale(:,e).*(flux1+flux2));
    
    flux1 = .5*diag([-1 1])*f3f(:,e);
    flux2 = sum((.5*diag([-1 1])*(-tf*VqPq)).*f3(rhoxf,rhoyf,uxf,uyf,Exf,Eyf),2);
    rhs3(:,e) = rx(:,e).*(Pq*sum(Dq.*f3(rhox,rhoy,ux,uy,Ex,Ey),2)) + Fscale(1,e)*(Pq*sum(Dfq.*f3(rhox,rhoy,ux,uy,Ex,Ey),2)) + LIFT*(Fscale(:,e).*(flux1+flux2));
    
end

global tau

% local lax penalty
cvel = sqrt(gamma*pM./rhoM);
lm   = (abs(uM) + cvel);
Lfc  = max(lm(mapP),lm);
% Lfc  = .5*(lm(mapP)+lm);

rhs1 = -2*rhs1 + .5*tau*LIFT*(Fscale.*Lfc.*(rhoP-rhoM));
rhs2 = -2*rhs2 + .5*tau*LIFT*(Fscale.*Lfc.*(mP-mM));
rhs3 = -2*rhs3 + .5*tau*LIFT*(Fscale.*Lfc.*(EP-EM));

end