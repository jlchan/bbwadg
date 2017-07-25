clear
Globals1D;

N = 4;
K1D = 8;
tau = 0;
FinalTime = 2;
CFL = .25;

r = JacobiGL(0,0,N);
% r = JacobiGQ(0,0,N);

[rq wq] = JacobiGL(0,0,N);
[rq wq] = JacobiGQ(0,0,N+1);

% % include boundary nodes
rq = [-1;rq;1];
wq = [0;wq;0];
Nq = length(rq);

% evaluate energy
rq2 = rq; wq2 = wq;
% [rq2 wq2] = JacobiGQ(0,0,N+1);
% [rq2 wq2] = JacobiGL(0,0,N);


[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

mapM = reshape(1:Nfp*Nfaces*K,Nfp*Nfaces,K);
mapP = mapM;
for e = 1:K
    mapP(:,e) = [(Nfp*Nfaces)*e-2, Nfp*Nfaces*(e+1)-1];  
end
% mapP(1) = 1; mapP(end) = Nfp*Nfaces*K;
mapP(1) = mapM(end); mapP(end) = mapM(1);
vmapP(1) = vmapM(end); vmapP(end) = vmapM(1);

V = Vandermonde1D(N,r);
Dr = GradVandermonde1D(N,r)/V;
Vq = Vandermonde1D(N,rq)/V;
Vf = Vandermonde1D(N,[-1 1])/V;
Vq2 = Vandermonde1D(N,rq2)/V;

M = Vq'*diag(wq)*Vq;
LIFT = M\Vf';

Vrq = GradVandermonde1D(N,rq)/V;
Pq = M\(Vq'*diag(wq));
Prq = M\(Vrq'*diag(wq));
xq =  Vandermonde1D(N,rq)/Vandermonde1D(N,JacobiGL(0,0,N))*x;
Pq(abs(Pq)<1e-8) = 0;

nxJ = nx.*Fscale;
du = @(uf) reshape(uf(mapP) - uf(mapM), Nfp*Nfaces, K);
Dh = @(u) rx.*(Dr*u) + .5*LIFT*(du(Vf*u).*nxJ);

D = zeros((N+1)*K);
uu = zeros(N+1,K);
for i = 1:(N+1)*K
    uu(i) = 1;
    tmp = Dh(uu);
    D(:,i) = tmp(:);
    uu(i) = 0;
end
D(abs(D)<1e-8) = 0;

% lift matrix
S = zeros((N+1)*K);
uu = zeros(N+1,K);
for i = 1:(N+1)*K
    uu(i) = 1;
    uuf = Vf*uu;
    du = uuf(mapP)-uuf(mapM);
    tmp = LIFT*reshape(du,Nfp*Nfaces,K);
    S(:,i) = tmp(:);
    uu(i) = 0;
end

rp = linspace(-1,1,100)';
Vp = Vandermonde1D(N,rp)/V;
xp = Vandermonde1D(N,rp)/Vandermonde1D(N,JacobiGL(0,0,N))*x;

Dx = kron(spdiag(rx(1,:)),Dr);
Dq = kron(speye(K),Vq)*D*kron(speye(K),Pq);
Dq(abs(Dq)<1e-8) = 0;
Dq = sparse(Dq);
PqG = kron(speye(K),Pq);
MM = kron(spdiag(J(1,:)),Vq'*diag(wq)*Vq);

VqG = kron(speye(K),Vq);
VpG = kron(speye(K),Vp);
xq = xq(:); xp = xp(:); 
wJq = diag(wq)*(Vq*J); wJq = wJq(:);
wJq2 = diag(wq2)*(Vq2*J); wJq2 = wJq2(:);
Vq2G = kron(speye(K),Vq2);

W = diag(wJq);

% chen/shu ops
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

% extract from all points
tM = zeros(Nfp*Nfaces*K,Nq*K);
tP = zeros(Nfp*Nfaces*K,Nq*K);
uu = zeros(Nq,K);
for i = 1:Nq*K
    uu(i) = 1;
    uf = uu([1; Nq],:);
    
    tmp = uf(mapM);
    tM(:,i) = tmp(:);
    tmp = uf(mapP);
    tP(:,i) = tmp(:);
    
    uu(i) = 0;
end

nx = nx(:);
Dh = .5*(VqG*Dx*PqG - VqG*(MM\(Dx'*VqG'*W))) + .5*VqG*L*(tM);
Fh = .5*tM'*diag(nx(:))*(tP - tM*VqG*PqG);
Pfh = MM\(VqG');

%% fluxes
g = 1;
opt = 1;
if opt==1 % shu's paper - equiv to gassner
    
    f1 = @(hL,hR,vL,vR) (hL.*vL + hR.*vR)/2;
    %f2 = @(hL,hR,vL,vR) .25*(hL.*vL.*vL + hR.*vR.*vL + hL.*vL.*vR) + .5*g.*(hL.*hR);
    f2 = @(hL,hR,vL,vR) .5*(hL.*vL + hR.*vR).*.5.*(vL+vR) + .5*g.*(hL.*hR);
    
elseif opt==2
    
    % % split form in gassner's paper?
    f1 = @(hL,hR,vL,vR) (hL.*vL + hR.*vR)/2;
    f2 = @(hL,hR,vL,vR) .5*( .5*(hL.*vL.*vL + hR.*vR.*vL + hL.*vL.*vR) + g*(hL.*hR));
    
elseif opt==3    
    
    % % EC flux function in Gassner's paper
    f1 = @(hL,hR,vL,vR) (hL + hR)/2.*(vL+vR)/3;
    f2 = @(hL,hR,vL,vR) .5*((.5*(vL+vR)).^2.*.5*(hL+hR) + .5*g.*(hL.^2 + hR.^2));   
    
end
%%
res1 = 0;
res2 = 0;

CN = (N+1)^2/2;
h = 2/K1D;
dt = CFL * h/CN;
Nsteps = ceil(FinalTime/dt); 
dt = FinalTime/Nsteps;

hex  = @(x) sin(pi*x);
hvex = @(x) cos(pi*x);
hex  = @(x) 2 + exp(-10*x.^2);
hvex = @(x) zeros(size(x));

% h = P*hex(xq);
h = hex(x(:));
% hv = P*hvex(xq);
hv = hvex(x(:));

figure
for i = 1:Nsteps
    for INTRK = 1:5
        
        % interpolate to quadrature
        hq = VqG*h;
        hvq = VqG*hv;        
        vq = hvq./hq;
        
        % project entropy variables
        q1 = VqG*PqG*(g*hq-vq.^2/2);
        q2 = VqG*PqG*(vq);
        
        % redefine flux variables
        hq = (q1+q2.^2/2)/g;
        vq = q2;
        
        [hx hy] = meshgrid(hq);
        [vx vy] = meshgrid(vq);
        
        %         d1 = PqG*(sum(Dq.*f1(hx,hy,vx,vy),2));
        %         rhs1 = 2*d1 + tau*S*h;
        %         d1 = PqG*(sum(Dq.*f2(hx,hy,vx,vy),2));
        %         rhs2 = 2*d1 + tau*S*hv;
        
        d1 = PqG*(sum(Dh.*f1(hx,hy,vx,vy),2)) + Pfh*(sum(Fh.*f1(hx,hy,vx,vy),2));
        rhs1 = 2*d1 + tau*S*h;
        d1 = PqG*(sum(Dh.*f2(hx,hy,vx,vy),2)) + Pfh*(sum(Fh.*f2(hx,hy,vx,vy),2));
        rhs2 = 2*d1 + tau*S*hv;
        
        % try to recover total energy = entropy
        if INTRK==5
            rhstest(i) = (g*hq - vq.^2/2)'*(wJq.*(VqG*rhs1)) + vq'*(wJq.*(VqG*rhs2));
        end
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        h = h + rk4b(INTRK)*res1;
        hv = hv + rk4b(INTRK)*res2;
        
    end;
    
        
    hq = Vq2G*h;
    hvq = Vq2G*hv;
    
%     vq = Vq2G*P*(hvq./hq);    
    vq = hvq./hq;
    
    energy(i) = wJq2(:)'*(hq.*vq.^2 + g*hq.^2)/2;
    
    if mod(i,25)==0
        hp = VpG*h;
        hvp = VpG*hv;        
        plot(xp,hp,'o')
        hold on
        plot(xp,hvp./hp,'x')
        title(sprintf('Time = %f',dt*i))
        axis([-1,1 -1 5])
        hold off
        drawnow
    end    
end

dS = max(energy)-min(energy)

figure(2)
semilogy(dt* (1:Nsteps),abs(energy-energy(1)),'--','linewidth',2)
hold on
title(sprintf('mean of rhs tested with h,v = %g, dS = %g\n',mean(rhstest),dS))

% axis([0 FinalTime energy(1)*.99 energy(1)*1.01])



