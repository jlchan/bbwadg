Globals1D
global rhoL rhoR mL mR EL ER gamma
global mapP mapB
global f1 f2 f3
global DNr Pq Lq
global Vf wq Vq

load Usnap_smooth.mat
% load Usnap_shock.mat

[U S V] = svd([reshape(Usnap,size(Usnap,1)/3,3*size(Usnap,2)) reshape(Vsnap,size(Vsnap,1)/3,3*size(Vsnap,2))],0);

s = diag(S);
Nmodes = 10;
tscale = 2;
global tau
tau = 1;

[Ur,Sr,~] = svd([U(:,1:Nmodes-1) ones(size(U,1),1)],0);
plot(x(:),Ur(:,end))
return

wJq = diag(wq)*(Vq*J);
wJq = wJq(:);

VqK = kron(eye(K),Vq)*Ur;
VNK = kron(eye(K),[Vq;Vf])*Ur;
MK = VqK'*diag(wJq)*VqK;
PK = MK\(VqK'*diag(wJq));

id = 1; %size(Usnap,2)-1;
U0 = reshape(Usnap(:,id),(N+1)*K,3);
V0 = reshape(Vsnap(:,id),(N+1)*K,3);
rho = VqK*PK*kron(eye(K),Vq)*U0(:,1);
m = VqK*PK*kron(eye(K),Vq)*U0(:,2);
E = VqK*PK*kron(eye(K),Vq)*U0(:,3);

v1 = V1(rho,m,E);
v2 = V2(rho,m,E);
v3 = V3(rho,m,E);

subplot(2,1,1)
hold on
plot(xq(:),rho,'--')
plot(xq(:),m,'--')
plot(xq(:),E,'--')

plot(x(:),U0,'o')

subplot(2,1,2)
hold on
plot(xq(:),v1,'-')
plot(xq(:),v2,'-')
plot(xq(:),v3,'-')
plot(x(:),V0,'o')

% figure
% semilogy(diag(S),'x')
% title(sprintf('energy in %d modes: %f\n',Nmodes,sum(s(1:Nmodes).^2)./sum(s(:).^2)))

%%
close all

% periodic BCs
mapP(1) = mapM(end); mapP(end) = mapM(1);
vmapP(1) = vmapM(end); vmapP(end) = vmapM(1);

dt = tscale*dt;
Nsteps = Nsteps/tscale;

rho = Pq*reshape(rho,Nq,K);
m = Pq*reshape(m,Nq,K);
E = Pq*reshape(E,Nq,K);

res1 = zeros(N+1,K);
res2 = zeros(N+1,K);
res3 = zeros(N+1,K);

for i = 1:Nsteps
    
    for INTRK = 1:5               
        
        rhoq = Vq*rho;
        mq = Vq*m;
        Eq = Vq*E;
        
        % interpolate to quadrature                
        if any(rhoq<1e-10)
            keyboard
        end
        
        % project entropy variables
        q1 = reshape(VNK*PK*V1(rhoq(:),mq(:),Eq(:)),(Nq+2),K); 
        q2 = reshape(VNK*PK*V2(rhoq(:),mq(:),Eq(:)),(Nq+2),K); 
        q3 = reshape(VNK*PK*V3(rhoq(:),mq(:),Eq(:)),(Nq+2),K); 
        
        % eval conservative vars and RHS
        rhoq = reshape(U1(q1,q2,q3),Nq+2,K);
        mq   = reshape(U2(q1,q2,q3),Nq+2,K);
        Eq   = reshape(U3(q1,q2,q3),Nq+2,K);
        [rhs1 rhs2 rhs3] = rhsEuler(rhoq,mq,Eq);        
        rhs1q = Vq*rhs1;
        rhs2q = Vq*rhs2;
        rhs3q = Vq*rhs3;
        rhs1 = Pq*reshape(VqK*PK*rhs1q(:),Nq,K);
        rhs2 = Pq*reshape(VqK*PK*rhs2q(:),Nq,K);
        rhs3 = Pq*reshape(VqK*PK*rhs3q(:),Nq,K);
        
        % project onto modes
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        rho = rho + rk4b(INTRK)*res1;
        m = m + rk4b(INTRK)*res2;
        E = E + rk4b(INTRK)*res3;
    end                
    
    if mod(i,5)==0 || i==Nsteps
        
        rhoq = Vq*rho;
        mq = Vq*m;
        Eq = Vq*E;
        
        plot(xq,rhoq,'b--','linewidth',2)
        hold on
        plot(xq,mq,'r--','linewidth',2)
        plot(xq,Eq,'k--','linewidth',2)
        
        title(sprintf('Time = %f, tstep %d out of %d',dt*i,i,Nsteps))

        hold off
        drawnow
        
    end        
end

%% final time

Uend = kron(eye(K),Vq)*reshape(Usnap(:,end-1),(N+1)*K,3);

hold on
% exact sol
plot(xq(:),Uend,'-','linewidth',2)

title(sprintf('Rel err = %g',norm(diag(wJq(:))*([rhoq(:) mq(:) Eq(:)] - Uend),'fro')/norm(diag(wJq(:))*Uend,'fro')))

%%

function [rhs1 rhs2 rhs3] = rhsEuler(rhoq,mq,Eq)

Globals1D
global rhoL rhoR mL mR EL ER gamma
global mapP mapB
global f1 f2 f3
global DNr Pq Lq
global Vf wq Vq

uq = mq./rhoq;

rhoM = rhoq(end-1:end,:);
uM = uq(end-1:end,:);
EM = Eq(end-1:end,:);
mM = rhoM.*uM;
rhoP = reshape(rhoM(mapP),Nfp*Nfaces,K);
uP = reshape(uM(mapP),Nfp*Nfaces,K);
EP = reshape(EM(mapP),Nfp*Nfaces,K);
mP = rhoP.*uP;
pM = (gamma-1)*(EM-.5*rhoM.*uM.^2);

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
    
    FS = f1(rhox,rhoy,ux,uy,Ex,Ey);
    FS = rx(1,e)*sum(DNr.*FS,2);
    rhs1(:,e) = [Pq Lq]*FS;
    
    FS = f2(rhox,rhoy,ux,uy,Ex,Ey);
    FS = rx(1,e)*sum(DNr.*FS,2);
    rhs2(:,e) = [Pq Lq]*FS;
    
    FS = f3(rhox,rhoy,ux,uy,Ex,Ey);
    FS = rx(1,e)*sum(DNr.*FS,2);
    rhs3(:,e) = [Pq Lq]*FS; 
end

global tau
% local lax penalty
cvel = sqrt(gamma*pM./rhoM);
lm   = (abs(uM) + cvel);
Lfc  = max(lm(mapP),lm);

d1 = Lfc.*((rhoP-rhoM));
d2 = Lfc.*((mP-mM));
d3 = Lfc.*((EP-EM));

f1f = -nx.*f1f + .5*tau*d1;
f2f = -nx.*f2f + .5*tau*d2;
f3f = -nx.*f3f + .5*tau*d3;
rhs1 = -2*rhs1 + Lq*(Fscale.*f1f);
rhs2 = -2*rhs2 + Lq*(Fscale.*f2f);
rhs3 = -2*rhs3 + Lq*(Fscale.*f3f);

end