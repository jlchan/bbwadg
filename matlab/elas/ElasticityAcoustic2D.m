function ElasticityAcoustic2D

% clear all, clear
clear -global *

Globals2D

K1D = 2;
N = 3;

FinalTime = .4;

global xp yp xq yq

% save acoustic mesh
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D,K1D/2);
VY = (VY+1)/2;
StartUp2D;

% cubature/plotting
Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
Vq = Vandermonde2D(N,rq,sq)/V;
Pq = V*V'*Vq'*diag(wq); % J's cancel out 

xpa = Vp*x; ypa = Vp*y;
Ka = K;
xq = Vq*x; yq = Vq*y;
meshAcoustic = getMesh();

% save elastic mesh
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D,K1D/2);
VY = (VY+1)/2-1;
StartUp2D;
xpe = Vp*x; ype = Vp*y; Ke = K;
xq = Vq*x; yq = Vq*y;
meshElastic = getMesh();

% unified plotting
xp(:,1:Ka) = xpa;
yp(:,1:Ka) = ypa;
xp(:,(1:Ke)+Ka) = xpe;
yp(:,(1:Ke)+Ka) = ype;

if 0
    % check
    getMesh(meshAcoustic);
    plot(x,y,'x')
    
    hold on
    getMesh(meshElastic);
    plot(x,y,'o')
    return
end

%% set up mesh coupling

global mapA vmapA mapE vmapE

getMesh(meshAcoustic);
vmapMA = vmapM;
mapTmp = reshape(mapM,Nfp*Nfaces,K);
xf = x(vmapM); yf = y(vmapM);
off = 1;
for e = 1:K
    mapf = reshape(mapTmp(:,e),Nfp,Nfaces);
    for f = 1:Nfaces
        if norm(yf(mapf(:,f)))<1e-8
            mapA(:,off) = mapf(:,f);
            off = off + 1;
        end
    end
end
xfa = xf(mapA); yfa = yf(mapA);

getMesh(meshElastic);
vmapME = vmapM;
mapTmp = reshape(mapM,Nfp*Nfaces,K);
xf = x(vmapM); yf = y(vmapM);
off = 1;
for e = 1:K
    mapf = reshape(mapTmp(:,e),Nfp,Nfaces);
    for f = 1:Nfaces
        if norm(yf(mapf(:,f)))<1e-8
            mapE(:,off) = mapf(:,f);
            off = off + 1;
        end
    end
end
xfe = xf(mapE); yfe = yf(mapE);

e = ones(1,Nfp);
for f = 1:size(xfe,2)
    D = (xfa(:,f)*e - (xfe(:,f)*e)').^2 + (yfa(:,f)*e - (yfe(:,f)*e)').^2;
    [idi idj] = find(D<1e-8);
    idM(:,f) = idi + (f-1)*Nfp;
    idP(:,f) = idj + (f-1)*Nfp;
end

mapA = mapA(idM);
mapE = mapE(idP);
vmapA = vmapMA(mapA);
vmapE = vmapME(mapE);

% getMesh(meshAcoustic);
% [~,p] = sort(x(vmapA(:)));
% getMesh(meshElastic);
% [~,p] = sort(x(vmapE(:)));
% vmapA = vmapA(p); vmapE = vmapE(p);
% mapA = mapA(p);   mapE = mapE(p);

if 0
    for i = 1:length(vmapB(:))
        getMesh(meshAcoustic);
        plot(x(vmapB(i)),y(vmapB(i)),'o')
        hold on
        axis([-1 1 -1 1])
        pause
    end
    return
end

if 0
    for i = 1:length(vmapA(:))
        getMesh(meshAcoustic);
        plot(x(vmapM(mapA(i))),y(vmapM(mapA(i))),'o')
        hold on
        getMesh(meshElastic);
        plot(x(vmapP(mapE(i))),y(vmapP(mapE(i))),'x')
        axis([-1 1 -1 1])
        pause
    end
end
% keyboard

%%

getMesh(meshElastic);
global Nfld mu lambda Vq Pq tau tauv taus useWADG
Nfld = 5; %(u1,u2,sxx,syy,sxy)

mu = 0;
lambda = 1;

useWADG = 0;
if useWADG
    mu = Vq*mu;
    lambda = Vq*lambda;
end

tau0 = 0;

tau = tau0;
tauv = tau0;
taus = tau0;

getMesh(meshAcoustic);
global cq
cq = ones(size(xq));

%% initial condition setup

% acoustic
getMesh(meshAcoustic);
x0 = 0;
y0 = .25;

pp = exp(-10^2*((x-x0).^2 + (y-y0).^2));
k = 1;
ww = (1-x).*(1+x).*y.*(1-y);%cos(k*pi*x/2).*cos(k*pi*y/2).*y;
z = zeros(Np, K);

Ua{1} = pp;
Ua{2} = z;
Ua{3} = z;

% elasticity
getMesh(meshElastic);
x0 = 0; y0 = -.25;
pp = exp(-10^2*((x-x0).^2 + (y-y0).^2));
z = zeros(Np, K);
Ue{1} = z;
Ue{2} = z;
Ue{3} = z;
Ue{4} = z;
Ue{5} = z;

%% compute eigenvalues
totaldofs = Np*Ka*3 + Np*Ke*5;

if 1 && totaldofs < 2000
    UA = zeros(Np*Ka,3);
    rhsa = zeros(Np*Ka,3);
    UE = zeros(Np*Ke,5);
    rhse = zeros(Np*Ke,5);
    
    A = zeros(totaldofs);
    
    % acoustic dofs
    for i = 1:Np*Ka*3
        UA(i) = 1;
        for fld = 1:3
            Ua{fld} = reshape(UA(:,fld),Np,Ka);
        end
        for fld = 1:5
            Ue{fld} = reshape(UE(:,fld),Np,Ke);
        end
        getMesh(meshAcoustic);
        rhsA = AcousRHS2D(Ua,Ue);
        getMesh(meshElastic);
        rhsE = ElasRHS2D(Ue,Ua);
        
        for fld = 1:3
            rhsa((1:Np*Ka) + (fld-1)*Np*Ka) = rhsA{fld};
        end
        for fld = 1:5
            rhse((1:Np*Ke) + (fld-1)*Np*Ke) = rhsE{fld};
        end
        A(:,i) = [rhsa(:); rhse(:)];
        
        UA(i) = 0;
        
        if mod(i,round(Np*Ka*3)/10)==0
            disp(sprintf('on acoustic dofs %d out of %d\n',i,Np*Ka*3));
        end
    end
    
    
    % elastic side
    for i = 1:Np*Ke*5
        UE(i) = 1;
        
        for fld = 1:3
            Ua{fld} = reshape(UA(:,fld),Np,Ka);
        end
        for fld = 1:5
            Ue{fld} = reshape(UE(:,fld),Np,Ke);
        end
        getMesh(meshAcoustic);
        rhsA = AcousRHS2D(Ua,Ue);
        getMesh(meshElastic);
        rhsE = ElasRHS2D(Ue,Ua);
        
        for fld = 1:3
            rhsa((1:Np*Ka) + (fld-1)*Np*Ka) = rhsA{fld};
        end
        for fld = 1:5
            rhse((1:Np*Ke) + (fld-1)*Np*Ke) = rhsE{fld};
        end
        A(:,i+Np*Ka*3) = [rhsa(:); rhse(:)];
        
        UE(i) = 0;
        
        if mod(i,round(Np*Ka*3)/10)==0
            disp(sprintf('on elastic dofs %d out of %d\n',i,Np*Ke*5));
        end
    end
    
    lam = eig(A);
    plot(lam,'o')
    title(sprintf('largest real part = %g\n',max(real(lam))))
    keyboard
end


%%

time = 0;
for fld = 1:3
    resA{fld} = zeros(Np,K);
end

for fld = 1:Nfld
    resE{fld} = zeros(Np,K);
end

% compute time step size
CN = (N+1)^2/2; % guessing...
dt = 2/(max(2*mu(:)+lambda(:))*CN*max(Fscale(:)));

% outer time step loop
tstep = 0;

M = inv(V*V');
wqJ = diag(wq)*(Vq*J);

while (time<FinalTime)
    
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        timeloc = time + rk4c(INTRK)*dt;
        
        getMesh(meshAcoustic);
        rhsA = AcousRHS2D(Ua,Ue);
        getMesh(meshElastic);
        rhsE = ElasRHS2D(Ue,Ua);
        
        % update both fields
        for fld = 1:3
            resA{fld} = rk4a(INTRK)*resA{fld} + dt*rhsA{fld};
            Ua{fld} = Ua{fld} + rk4b(INTRK)*resA{fld};
        end
        
        for fld = 1:Nfld
            resE{fld} = rk4a(INTRK)*resE{fld} + dt*rhsE{fld};
            Ue{fld} = Ue{fld} + rk4b(INTRK)*resE{fld};
        end
        
        pp(:,1:Ka) = Ua{1};
        pp(:,(1:Ke)+Ka) = (Ue{3}+Ue{4})/2;
% pp(:,1:Ka) = Ua{3};
% pp(:,(1:Ke)+Ka) = Ue{2};

    end
    
    if 1 && mod(tstep,10)==0
        clf
        vv = Vp*pp;
        color_line3(xp,yp,vv,vv,'.');
        axis tight
        
        %         Ue = U;
        %         Snx = Ue{3}(vmapE).*nx(mapA) + Ue{4}(vmapE).*ny(mapA);
        %         Sny = Ue{4}(vmapE).*nx(mapA) + Ue{5}(vmapE).*ny(mapA);
        %         nSn = Snx.*nx(mapA) + Sny.*ny(mapA);
        %         plot(x(vmapA),nSn,'o')
        %         axis([-1 1 -1 1])
        
        title(sprintf('time = %f',time));
        colorbar;
        drawnow
    end
    
    time = time+dt;
    tstep = tstep+1;
    
end

clf
vv = Vp*pp;
color_line3(xp,yp,vv,vv,'.');
axis tight
title(sprintf('time = %f',time));
colorbar;

keyboard
return



function [rhs] = ElasRHS2D(U,Ua)

Globals2D;

global Nfld mu lambda Vq Pq tauv taus useWADG
global C11 C12 C13 C22 C23 C33
global mapA vmapA mapE vmapE

% Define field differences at faces
for fld = 1:Nfld
    u = U{fld};
    
    % compute jumps
    dU{fld} = zeros(Nfp*Nfaces,K);
    dU{fld}(:) = u(vmapP)-u(vmapM);
    
    ur = Dr*u;
    us = Ds*u;
    Ux{fld} = rx.*ur + sx.*us;
    Uy{fld} = ry.*ur + sy.*us;
end

divSx = Ux{3} + Uy{5}; % d(Sxx)dx + d(Sxy)dy
divSy = Ux{5} + Uy{4}; % d(Sxy)dx + d(Syy)dy
du1dx = Ux{1}; % du1dx
du2dy = Uy{2}; % du2dy
du12dxy = Ux{2} + Uy{1}; % du2dx + du1dy

% velocity fluxes
Snx = nx.*dU{3} + ny.*dU{5};
Sny = nx.*dU{5} + ny.*dU{4};

opt = 1;
if opt==1 % traction BCs
    Snx(mapB) = -2*(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    Sny(mapB) = -2*(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
elseif opt==2 % basic ABCs
    Snx(mapB) = -(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    Sny(mapB) = -(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
    dU{1}(mapB) = -U{1}(vmapB);
    dU{2}(mapB) = -U{2}(vmapB);
elseif opt==3 % zero velocity
    dU{1}(mapB) = -2*U{1}(vmapB);
    dU{2}(mapB) = -2*U{2}(vmapB);
    
end

% stress fluxes
nUx = dU{1}.*nx;
nUy = dU{2}.*ny;
nUxy = dU{2}.*nx + dU{1}.*ny;

% evaluate central fluxes
fc{1} = Snx;
fc{2} = Sny;
fc{3} = nUx;
fc{4} = nUy;
fc{5} = nUxy;

% penalization terms - reapply An
fp{1} = nx.*fc{3} + ny.*fc{5};
fp{2} = nx.*fc{5} + ny.*fc{4};
fp{3} = fc{1}.*nx;
fp{4} = fc{2}.*ny;
fp{5} = fc{2}.*nx + fc{1}.*ny;

flux = cell(5,1);
for fld = 1:Nfld
    flux{fld} = zeros(Nfp*Nfaces,K);
    if fld < 3
        taufld = tauv;
    else
        taufld = taus;
    end
    flux{fld}(:) = fc{fld}(:) + taufld.*fp{fld}(:);
end

% compute acoustic-elastic interface flux (central)
nxf = nx(mapE); nyf = ny(mapE);
dU1 = (Ua{2}(vmapA) - U{1}(vmapE));
dU2 = (Ua{3}(vmapA) - U{2}(vmapE));
dUx = dU1.*nxf;
dUy = dU2.*nyf;
dUxy = dU2.*nxf + dU1.*nyf;

pf = Ua{1}(vmapA);
Snx = nxf.*U{3}(vmapE) + nyf.*U{5}(vmapE);
Sny = nxf.*U{5}(vmapE) + nyf.*U{4}(vmapE);
dSx = (pf.*nxf - Snx); 
dSy = (pf.*nyf - Sny); 

% fc{1} = Snx + (nx.*fc{3} + ny.*fc{5});
% fc{2} = Sny + (nx.*fc{5} + ny.*fc{4});
% fc{3} = nUx + (fc{1}.*nx);
% fc{4} = nUy + (fc{2}.*ny);
% fc{5} = nUxy + (fc{2}.*nx + fc{1}.*ny);

% flux{1}(mapE) = dSx + tauv.*(dUx.*nxf + dUxy.*nyf);
% flux{2}(mapE) = dSy + tauv.*(dUxy.*nxf + dUy.*nyf);
% flux{3}(mapE) = dUx + taus.*(dSx.*nxf);
% flux{4}(mapE) = dUy + taus.*(dSy.*nyf);
% % flux{5}(mapE) = dUx.*nyf + dUy.*nxf + taus.*(dSx.*nyf + dSy.*nxf);
% flux{5}(mapE) = dUxy + taus.*(dSx.*nyf + dSy.*nxf);

% compute right hand sides of the PDE's
rr{1} =  divSx   +  LIFT*(.5*Fscale.*flux{1});
rr{2} =  divSy   +  LIFT*(.5*Fscale.*flux{2});
rr{3} =  du1dx   +  LIFT*(.5*Fscale.*flux{3});
rr{4} =  du2dy   +  LIFT*(.5*Fscale.*flux{4});
rr{5} =  du12dxy +  LIFT*(.5*Fscale.*flux{5});

if useWADG
    for fld = 3:Nfld
        rr{fld} = Vq*rr{fld};
    end
    rhs{1} = rr{1};
    rhs{2} = rr{2};
    rhs{3} = Pq*((2*mu+lambda).*rr{3} + lambda.*rr{4});
    rhs{4} = Pq*(lambda.*rr{3} + (2*mu+lambda).*rr{4});
    rhs{5} = Pq*((mu) .* rr{5});
else
    rhs{1} = rr{1};
    rhs{2} = rr{2};
    rhs{3} = ((2*mu+lambda).*rr{3} + lambda.*rr{4});
    rhs{4} = (lambda.*rr{3} + (2*mu+lambda).*rr{4});
    rhs{5} = (mu) .* rr{5};
end

% for fld = 1:5
%     rhs{fld} = 0*rhs{fld};
% end


return;



function rhs = AcousRHS2D(U,Ue)

Globals2D;
global cq Vq Pq
global mapA vmapA mapE vmapE

p = U{1}; u = U{2}; v = U{3};

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

ndotdU = nx.*du + ny.*dv;

% Impose reflective boundary conditions (p+ = -p-)
dp(mapB) = -2*p(vmapB);
ndotdU(mapB) = 0;

global tau
fluxp = ndotdU + tau*dp;
fluxu = (dp + tau*ndotdU).*nx;
fluxv = (dp + tau*ndotdU).*ny;

% compute acoustic-elastic interface flux
nxf = nx(mapA); nyf = ny(mapA);
pf = p(vmapA); 
uf = u(vmapA); vf = v(vmapA);

dUn = (Ue{1}(vmapE)-uf).*nxf + (Ue{2}(vmapE)-vf).*nyf; 
Snx = Ue{3}(vmapE).*nxf + Ue{5}(vmapE).*nyf;
Sny = Ue{5}(vmapE).*nxf + Ue{4}(vmapE).*nyf;

% fluxp(mapA) = dUn + tau*((Snx.*nxf + Sny.*nyf)-pf);
% fluxu(mapA) = ((Snx - pf.*nxf) + tau*dUn).*nxf;
% fluxv(mapA) = ((Sny - pf.*nyf) + tau*dUn).*nyf;

pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
rhs{1} =  divU + LIFT*(.5*Fscale.*fluxp);
rhs{2} =  dpdx + LIFT*(.5*Fscale.*fluxu);
rhs{3} =  dpdy + LIFT*(.5*Fscale.*fluxv);

% rhs{1} = Pq*(cq.*(Vq*rhs{1}));
% for fld = 1:3
%     rhs{fld} = 0*rhs{fld};
% end


return

