function ElasticityAcoustic2D

% clear all, clear
clear -global *

Globals2D

K1D = 16;
N = 5;
c_flag = 0;
FinalTime = 2;

global xp yp xq yq

% save acoustic mesh
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D,K1D/2);
VY = (VY+1)/2;
StartUp2D;

% cubature/plotting
[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xpa = Vp*x; ypa = Vp*y; 
Ka = K;

Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
Pq = V*V'*Vq'*diag(wq); % J's cancel out
Mref = inv(V*V');

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

%%

getMesh(meshElastic);
global Nfld mu lambda Vq Pq tauv taus useWADG
Nfld = 5; %(u1,u2,sxx,syy,sxy)

mu = ones(size(x));
lambda = ones(size(x));

useWADG = 0;
if useWADG
    mu = Vq*mu;
    lambda = Vq*lambda;
end

tau0 = 1;
tauv = tau0*ones(size(x));
taus = tau0*ones(size(x));

getMesh(meshAcoustic);
global cq
cq = ones(size(xq));

%% initial condition setup

% acoustic
getMesh(meshAcoustic);
x0 = 0; y0 = .5;
pp = exp(-25^2*((x-x0).^2 + (y-y0).^2));
z = zeros(Np, K);
Ua{1} = 0*pp;
Ua{2} = z;
Ua{3} = z;

% elasticity
getMesh(meshElastic);
x0 = 0; y0 = -.5;
pp = exp(-25^2*((x-x0).^2 + (y-y0).^2));
z = zeros(Np, K);
Ue{1} = z;
Ue{2} = pp;
Ue{3} = z;
Ue{4} = z;
Ue{5} = z;

%%
getMesh(meshAcoustic);
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
        pp(:,1:Ka) = Ua{1};
        
        for fld = 1:Nfld
            resE{fld} = rk4a(INTRK)*resE{fld} + dt*rhsE{fld};
            Ue{fld} = Ue{fld} + rk4b(INTRK)*resE{fld};
        end
        pp(:,(1:Ke)+Ka) = Ue{3}+Ue{4};
        
    end
    
    
%     if abs(time-.05)<2*dt
%         getMesh(meshAcoustic);
%         p = Ua{1}; u = Ua{2}; v = Ua{3};
%         nxf = nx(mapA);
%         nyf = ny(mapA);
%         Snx = Ue{3}(vmapE).*nxf + Ue{5}(vmapE).*nyf;
%         Sny = Ue{5}(vmapE).*nxf + Ue{4}(vmapE).*nyf;
%         % nSn = Snx.*nxf + Sny.*nyf;
%         dUn = (Ue{1}(vmapE)-u(vmapA)).*nxf + (Ue{2}(vmapE)-v(vmapA)).*nyf;
%         fluxp(mapA) = dUn;
%         fluxu(mapA) = (Snx - p(vmapA).*nxf) + 2*p(vmapA).*nxf;
%         fluxv(mapA) = (Sny - p(vmapA).*nyf) + 2*p(vmapA).*nyf ;        
%         aflux = fluxp(mapA).*p(vmapA) + ...
%             fluxu(mapA).*u(vmapA) + ...
%             fluxv(mapA).*v(vmapA);
%         
%         v2f = Ue{2}(vmapE);
%         pf = p(vmapA);
%         uf = u(vmapA);
%         vf = v(vmapA);
%         aflux2 = (v2f-vf).*nyf.*pf + (Snx + pf.*nxf).*uf + (Sny + pf.*nyf).*vf;
% %         keyboard
% 
%         getMesh(meshElastic);
%         nxf = nx(mapE);
%         nyf = ny(mapE);
%         dUx = (Ua{2}(vmapA) - Ue{1}(vmapE));
%         dUy = (Ua{3}(vmapA) - Ue{2}(vmapE));
%         dUn = dUx.*nxf + dUy.*nyf;
%         pf = Ua{1}(vmapA);
%         nSx = nxf.*Ue{3}(vmapE) + nyf.*Ue{5}(vmapE);
%         nSy = nxf.*Ue{5}(vmapE) + nyf.*Ue{4}(vmapE);
%         
%         flux{1}(mapE) = (pf.*nxf - nSx) + 2*nSx;
%         flux{2}(mapE) = (pf.*nyf - nSy) + 2*nSy;
%         flux{3}(mapE) = dUx.*nxf;
%         flux{4}(mapE) = dUy.*nyf;
%         flux{5}(mapE) = dUx.*nyf + dUy.*nxf;
%         
%         eflux = flux{1}(mapE).*Ue{1}(vmapE) + ...
%             flux{2}(mapE).*Ue{2}(vmapE) + ...
%             flux{3}(mapE).*Ue{3}(vmapE) + ...
%             flux{4}(mapE).*Ue{4}(vmapE) + ...
%             flux{5}(mapE).*Ue{5}(vmapE);
%         
%         v1f = Ue{1}(vmapE);
%         v2f = Ue{2}(vmapE);
%         sxxf = Ue{3}(vmapE);
%         syyf = Ue{4}(vmapE);
%         sxyf = Ue{5}(vmapE);
%         
%         %aflux2 = (v2f-vf).*nyf.*pf + (Snx + pf.*nxf).*uf + (Sny+pf.*nyf).*vf;
%         eflux2 = (pf.*nxf + nSx).*v1f + (pf.*nyf + nSy).*v2f + ...
%             + dUx.*nxf.*sxxf + dUy.*nyf.*syyf + dUx.*nyf.*sxyf        
%         
%         keyboard
%     end
    
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
nSx = nx.*dU{3} + ny.*dU{5};
nSy = nx.*dU{5} + ny.*dU{4};

opt=1;
if opt==1 % traction BCs
    nSx(mapB) = -2*(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -2*(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
elseif opt==2 % basic ABCs
    nSx(mapB) = -(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
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
fc{1} = nSx;
fc{2} = nSy;
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
        taufld = tauv(vmapM);
    else
        taufld = taus(vmapM);
    end
    flux{fld}(:) = fc{fld}(:) + taufld.*fp{fld}(:);
end

% compute acoustic-elastic interface flux (central)
nxf = nx(mapE);
nyf = ny(mapE);
dUx = (Ua{2}(vmapA) - U{1}(vmapE));
dUy = (Ua{3}(vmapA) - U{2}(vmapE));
dUn = dUx.*nxf + dUy.*nyf;
pf = Ua{1}(vmapA);
nSx = nxf.*U{3}(vmapE) + nyf.*U{5}(vmapE);
nSy = nxf.*U{5}(vmapE) + nyf.*U{4}(vmapE);
dSx = (pf.*nxf - nSx);
dSy = (pf.*nyf - nSy);

% flux{1}(mapE) = dSx;
% flux{2}(mapE) = dSy;
% flux{3}(mapE) = dUx.*nxf;
% flux{4}(mapE) = dUy.*nyf;
% flux{5}(mapE) = dUx.*nyf + dUy.*nxf;

flux{1}(mapE) = dSx + tauv(mapE).*dUn.*nxf;
flux{2}(mapE) = dSy + tauv(mapE).*dUn.*nyf;
flux{3}(mapE) = dUx.*nxf + taus(mapE).*(dSx.*nxf);
flux{4}(mapE) = dUy.*nyf + taus(mapE).*(dSy.*nyf);
flux{5}(mapE) = dUx.*nyf + dUy.*nxf + taus(mapE).*(dSx.*nyf+dSy.*nxf);

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

tau = 1;
fluxp = tau*dp + ndotdU;
fluxu = (tau*ndotdU + dp).*nx;
fluxv = (tau*ndotdU + dp).*ny;

% compute acoustic-elastic interface flux
nxf = nx(mapA);
nyf = ny(mapA);
Snx = Ue{3}(vmapE).*nxf + Ue{5}(vmapE).*nyf;
Sny = Ue{5}(vmapE).*nxf + Ue{4}(vmapE).*nyf;
dUn = (Ue{1}(vmapE)-u(vmapA)).*nxf + (Ue{2}(vmapE)-v(vmapA)).*nyf;
% fluxp(mapA) = dUn;
% fluxu(mapA) = (Snx - p(vmapA).*nxf);
% fluxv(mapA) = (Sny - p(vmapA).*nyf);
fluxp(mapA) = tau*((Snx.*nxf + Sny.*nyf) - p(vmapA)) + dUn;
fluxu(mapA) = tau*dUn.*nx(mapA) + (Snx - p(vmapA).*nx(mapA));
fluxv(mapA) = tau*dUn.*ny(mapA) + (Sny - p(vmapA).*ny(mapA));

pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
rhs{1} =  divU + LIFT*(.5*Fscale.*fluxp);
rhs{2} =  dpdx + LIFT*(.5*Fscale.*fluxu);
rhs{3} =  dpdy + LIFT*(.5*Fscale.*fluxv);

rhs{1} = Pq*(cq.*(Vq*rhs{1}));


return

