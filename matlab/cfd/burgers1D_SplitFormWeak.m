clear
Globals1D;

N = 4;
K1D = 32;
CFL = .125;
T = 2;

[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
StartUp1D;
vmapP(1) = vmapM(end); vmapP(end) = vmapM(1);
mapM = reshape(1:Nfp*Nfaces*K,Nfp*Nfaces,K);
mapP = mapM;
for e = 1:K
    mapP(:,e) = [(Nfp*Nfaces)*e-2, Nfp*Nfaces*(e+1)-1];
end
mapM = mapM(:); mapP = mapP(:);
mapP(1) = mapM(end); mapP(end) = mapM(1);


r = JacobiGL(0,0,N);
% r = JacobiGQ(0,0,N);
[rq wq] = JacobiGL(0,0,N);
[rq wq] = JacobiGQ(0,0,N+1);

rq = [-1;rq;1];
wq = [0;wq;0];


V = Vandermonde1D(N,r);
Dr = GradVandermonde1D(N,r)/V;
Vq = Vandermonde1D(N,rq)/V;
Vf = Vandermonde1D(N,[-1 1])/V;

M = Vq'*diag(wq)*Vq;
LIFT = M\Vf';

Pq = M\(Vq'*diag(wq));
Prq = M\((Vq*Dr)'*diag(wq));
xq =  Vandermonde1D(N,rq)/Vandermonde1D(N,JacobiGL(0,0,N))*x;

rp = linspace(-1,1,100)';
Vp = Vandermonde1D(N,rp)/V;
xp = Vandermonde1D(N,rp)/Vandermonde1D(N,JacobiGL(0,0,N))*x;

x = x(:); xq = xq(:); xp = xp(:);
wJq = diag(wq)*(Vq*J); wJq = wJq(:);

Dx = kron(spdiag(rx(1,:)),Dr);
DTx = kron(spdiag(rx(1,:)),M\(Dr'*M));
PqG = kron(speye(K),Pq);
PrqG = kron(spdiag(rx(1,:)),Prq);
VqG = kron(speye(K),Vq);
VpG = kron(speye(K),Vp);
VfG = kron(speye(K),Vf);
MM = kron(spdiag(J(1,:)),Vq'*diag(wq)*Vq);

Dq = VqG*Dx*PqG;
DTq = VqG*PrqG;

% global lift
L = zeros((N+1)*K,Nfp*Nfaces*K);
uu = zeros(Nfp*Nfaces,K);
for i = 1:Nfp*Nfaces*K
    uu(i) = 1;
    tmp = LIFT*(reshape(uu,Nfp*Nfaces,K).*nx.*Fscale);
    L(:,i) = tmp(:);
    uu(i) = 0;
end
L(abs(L)<1e-8) = 0; L = sparse(L);

tM = zeros(Nfp*Nfaces*K);
tP = zeros(Nfp*Nfaces*K);
uu = zeros(Nfp*Nfaces,K);
for i = 1:Nfp*Nfaces*K
    uu(i) = 1;    
    
    tmp = uu(mapM);
    tM(:,i) = tmp(:);        
    tmp = uu(mapP);
    tP(:,i) = tmp(:);
    
    uu(i) = 0;
end

Nq = length(rq);
tq = zeros(Nfp*Nfaces*K,Nq*K);
uu = zeros(Nq,K);
for i = 1:Nq*K
    uu(i) = 1;    
    
    tmp = uu([1 Nq],:);
    tq(:,i) = tmp(:);
    
    uu(i) = 0;
end

% global DG deriv
Dhq = Dq + VqG*L*.5*(tP-tM)*tq*VqG*PqG;
Dhq2 = (Dq-DTq) + VqG*L*tP*tq*VqG*PqG;
%Dwq = .5*((Dq-DTq) + VqG*L*((tP+tM)*tq - tM*tq*VqG*PqG));
Dwq = .5*((Dq-DTq) + VqG*L*(tP*tq));

% Dw = Dq + VqG*(L*(.5*(tP+tM) - tM*VqG*PqG));
% Ds = Dq + VqG*L*.5*(tP-tM);

W = diag(wJq(:));

% norm(W*Dw-Ds'*W,'fro')
% return

fS = @(ux,uy) (ux.^2 + ux.*uy + uy.^2)/6;
%f = @(ux,uy) (ux.^2 + ux.*uy)/6;
f = @(ux,uy) (ux.^2 + uy.^2)/2;


resu = 0;

CN = (N+1)^2/2;
h = 2/K1D;
dt = CFL * h/CN;
Nsteps = ceil(T / dt);
dt = T/Nsteps;

u = PqG*sin(pi*xq);

% u = rand(Np*K,1);
% uM = VfG*u; uP = uM(mapP);
% return

% uq = VqG*u;
% [ux uy] = meshgrid(uq(:));
% PqG*sum((VqG*L*.5*(tP)*PqG).*f(ux,uy),2)
% return

figure
for i = 1:Nsteps
    for INTRK = 1:5
        uq = VqG*u;
        [ux uy] = meshgrid(uq(:));
        uM = VfG*u; uP = uM(mapP);
                        
%         rhsu = PqG*(sum((2*Dq).*f(ux,uy),2)) + L*fS(uM,uP) - (L*(VfG*PqG*uq.^2 + uM.^2/2))/3;        
%         rhsu = (-DTx*PqG*uq.^2 + PqG*(uq.*(Dq*uq)))/3 + L*(fS(uM,uP) - uM.^2/6);
        %rhsu = PqG*(DTq*uq.^2 + uq.*(Dq*uq))/3 + L*(uM.^2 + uM.*uP + uP.^2 - uM.^2)/6;
%         rhsu = PqG*(DTq*uq.^2 + uq.*(DTq*uq) + Dq*(uq.^2) + uq.*(Dq*uq))/6 + L*(uM.^2 + uM.*uP + uP.^2)/6 + ...
%             PqG*(uq.*(VqG*L*uM) - VqG*L*(uM.^2 + VfG*PqG*uq.^2))/6;
        
        %rhsu = -(DTx*PqG*uq.^2 + PqG*(uq.*(VqG*PrqG*uq)))/3 + L*(fS(uM,uP) - uM.^2/6) + PqG*(uq.*(VqG*L*uM/3));
        rhsu = -(DTx*PqG*uq.^2 + PqG*(uq.*(VqG*PrqG*uq)))/3 + L*((uM+uP).*uP/6) + PqG*(uq.*(VqG*L*uM/3));
                        
%         rhsu = PqG*(sum((Dq+DTq).*f(ux,uy),2)) + L*(fS(uM,uP) - uM.^2/2 + uM.^2/6 + uM.^2/6); % accurate for SEM only        
        
       
%         rhsu = 2*PqG*(sum(Dhq.*f(ux,uy),2));     
%         rhsu = 2*PqG*(sum(Dwq.*f(ux,uy),2));
  
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u + rk4b(INTRK)*resu;
    end;
    
    
    energy(i) = sum(sum(wJq.*(VqG*u).^2));
    if mod(i,10)==0
        plot(xp,VpG*u,'o')
        axis([-1,1 -4 4])
        drawnow
    end
end
% return
figure(2)
semilogy(dt*(1:Nsteps),energy,'o--')

dS = max(energy)-min(energy)
hold on
title(sprintf('dS = %g\n',dS))
