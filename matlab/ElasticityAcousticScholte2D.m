function ElasticityAcoustic2D

% clear all, clear
clear -global *

Globals2D

K1D = 4;
N = 4;
FinalTime = .5;

global xp yp xq yq Vq Pq

a = 2; % stretch factor

% save acoustic mesh
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D,K1D/2*a);
VY = (VY+1)/2;
VY = VY*a;
StartUp2D;

% BuildPeriodicMaps2D(2,-1); % in x-direction only

% cubature/plotting
[rp sp] = EquiNodes2D(100); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xpa = Vp*x; ypa = Vp*y;
Ka = K;

Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
Pq = V*V'*Vq'*diag(wq); % J's cancel out
Mref = inv(V*V');

xq = Vq*x; yq = Vq*y; wJq = diag(wq)*(Vq*J);
meshAcoustic = getMesh();

% save elastic mesh
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D,K1D/2*a);
VY = (VY+1)/2-1;
VY = VY*a;
StartUp2D;
% BuildPeriodicMaps2D(2,0); % in x-direction only

xpe = Vp*x; ype = Vp*y; Ke = K;
xq = Vq*x; yq = Vq*y; wJq = diag(wq)*(Vq*J);
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
    %     plot(x(vmapB),y(vmapB),'x')
    
    hold on
    getMesh(meshElastic);
    %     plot(x,y,'o')
    %     plot(x(vmapB),y(vmapB),'o')
    return
end

%% set up domain coupling connectivity

global mapA vmapA mapE vmapE

mapA = zeros(N+1,K1D);
mapE = zeros(N+1,K1D);
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

% match faces
e = ones(1,Nfp);
for f = 1:size(xfe,2)
    D = (xfa(:,f)*e - (xfe(:,f)*e)').^2 + (yfa(:,f)*e - (yfe(:,f)*e)').^2;
    [idi, idj] = find(D<1e-8);
    idM(:,f) = idi + (f-1)*Nfp;
    idP(:,f) = idj + (f-1)*Nfp;
end

mapA = mapA(idM);
mapE = mapE(idP);
vmapA = vmapMA(mapA);
vmapE = vmapME(mapE);
% keyboard

if 0
    % check
    getMesh(meshAcoustic);
    xf = x(vmapM);
    yf = y(vmapM);
    plot(xf(mapA),yf(mapA),'x')
    %     plot(x(vmapA),y(vmapA),'x')
    %     plot(x(vmapB),y(vmapB),'x')
    
    hold on
    getMesh(meshElastic);
    %     plot(x,y,'o')
    plot(x(vmapE),y(vmapE),'o')
    return
end


%% set up material parameters

getMesh(meshElastic);
global Nfld mu lambda tauv taus useWADG
Nfld = 5; % (u1,u2,sxx,syy,sxy)

lambda = 1;
mu = 1;

getMesh(meshAcoustic);
global cq
cq = ones(size(x)); % mu = 0 here
% cq = sqrt(2*mu+lambda); % cp

global tau
tau = 1;
tauv = tau*ones(size(x));
taus = tau*ones(size(x));


%% initial condition setup

% % acoustic
% getMesh(meshAcoustic);
% x0 = 0; y0 = .5;
% pp = exp(-10^2*((x-x0).^2 + (y-y0).^2));
% z = zeros(Np, K);
% Ua{1} = pp;
% Ua{2} = z;
% Ua{3} = z;
% 
% % elasticity
% getMesh(meshElastic);
% x0 = 0; y0 = -.5;
% pp = exp(-10^2*((x-x0).^2 + (y-y0).^2));
% z = zeros(Np, K);
% Ue{1} = z;
% Ue{2} = z;
% Ue{3} = z;
% Ue{4} = z;
% Ue{5} = z;

%% exact sol scholte

if 1
    mu1 = 0; % acoustic
    mu2 = mu; % elastic
    
    c1p = sqrt(2*mu1+lambda); % rho = 1
    c2p = sqrt(2*mu2+lambda);
    c2s = sqrt(mu2);
    
    c = 0.7110017230197;
    w = 2;
    k = w/c;
    B1 =  -1i*0.3594499773037;
    B2 =  -1i*0.8194642725978;
    B3 = 1;
    b1p = sqrt(1-c^2/c1p^2);
    b2p = sqrt(1-c^2/c2p^2);
    b2s = sqrt(1-c^2/c2s^2);
    
    global v1a v2a v1b v2b 
    v1a = @(x,y,t) real(B1.*k.*w.*exp(-b1p.*k.*y).*exp(k.*x.*1i - t.*w.*1i));
    v2a = @(x,y,t) real(B1.*b1p.*k.*w.*exp(-b1p.*k.*y).*exp(k.*x.*1i - t.*w.*1i).*1i);
    v1b = @(x,y,t) real(-k.*w.*exp(k.*x.*1i - t.*w.*1i).*(B2.*exp(b2p.*k.*y).*1i - B3.*b2s.*exp(b2s.*k.*y)).*1i);
    v2b = @(x,y,t) real(-k.*w.*exp(k.*x.*1i - t.*w.*1i).*(B2.*b2p.*exp(b2p.*k.*y) + B3.*exp(b2s.*k.*y).*1i).*1i);
    u1ax = @(x,y,t) real(-B1.*k^2.*exp(-b1p.*k.*y).*exp(k.*x.*1i - t.*w.*1i));
    u2ay = @(x,y,t) real(B1.*b1p^2.*k^2.*exp(-b1p.*k.*y).*exp(k.*x.*1i - t.*w.*1i));
    u12axy = @(x,y,t) real(-B1.*b1p.*k^2.*exp(-b1p.*k.*y).*exp(k.*x.*1i - t.*w.*1i).*2i);
    u1bx = @(x,y,t) real(k^2.*exp(k.*x.*1i - t.*w.*1i).*(B2.*exp(b2p.*k.*y).*1i - B3.*b2s.*exp(b2s.*k.*y)).*1i);
    u2by = @(x,y,t) real(k^2.*exp(k.*x.*1i - t.*w.*1i).*(B2.*b2p^2.*exp(b2p.*k.*y) + B3.*b2s.*exp(b2s.*k.*y).*1i));
    u12bxy = @(x,y,t) real(-exp(k.*x.*1i - t.*w.*1i).*(B3.*b2s^2.*k^2.*exp(b2s.*k.*y) - B2.*b2p.*k^2.*exp(b2p.*k.*y).*1i) + k.*exp(k.*x.*1i - t.*w.*1i).*(B2.*b2p.*k.*exp(b2p.*k.*y) + B3.*k.*exp(b2s.*k.*y).*1i).*1i);
    
    getMesh(meshAcoustic);
    Ua{1} = Pq*(u1ax(xq,yq,0) + u2ay(xq,yq,0)); % lambda*tr(E) = p?
    Ua{2} = Pq*v1a(xq,yq,0);
    Ua{3} = Pq*v2a(xq,yq,0);
    
    getMesh(meshElastic);
    Ue{1} = Pq*v1b(xq,yq,0);
    Ue{2} = Pq*v2b(xq,yq,0);
    Ue{3} = Pq*((2*mu2+lambda) .* u1bx(xq,yq,0) + lambda.*u2by(xq,yq,0));
    Ue{4} = Pq*(lambda.*u1bx(xq,yq,0) + (2*mu2+lambda) .* u2by(xq,yq,0));
    Ue{5} = Pq*(mu2 .* u12bxy(xq,yq,0));
    
    if 0
        figure
        pp(:,1:Ka) = Ua{1};
        pp(:,(1:Ke)+Ka) = Ue{3};
        vv = Vp*pp;
        color_line3(xp,yp,vv,vv,'.');
        colorbar
        %             return
        
        % check flux condition
        t = 0;
        nxx = 0; nyy = 1;
        xx = -1:.1:1;
        Sxx = (2*mu2+lambda) .* u1bx(xx,0,t) + lambda.*u2by(xx,0,t);
        Syy = lambda.*u1bx(xx,0,t) + (2*mu2+lambda) .* u2by(xx,0,t);
        Sxy = mu2 * u12bxy(xx,0,t);
        Snx = nxx*Sxx + nyy*Sxy;
        Sny = nxx*Sxy + nyy*Syy;
        v1 = v1b(xx,0,t);
        v2 = v2b(xx,0,t);
        p = u1ax(xx,0,t) + u2ay(xx,0,t);
        u = v1a(xx,0,t);
        v = v2a(xx,0,t);
        norm(Snx - p*nxx)
        norm(Sny - p*nyy)
        norm((v1 - u)*nxx)
        norm((v2 - v)*nyy)
        
        keyboard
        
        
        for t = 0:.1:.25
            vv(:,1:Ka) = u1ax(xpa,ypa,t) + u2ay(xpa,ypa,t);
            vv(:,(1:Ke)+Ka) = (2*mu2+lambda) .* u1bx(xpe,ype,t) + lambda.*u2by(xpe,ype,t);
            clf
            color_line3(xp,yp,vv,vv,'.');
            colorbar
            drawnow
        end
        return
    end
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
dt = .5/(CN*max(Fscale(:)));

% outer time step loop
tstep = 0;

M = inv(V*V');
wqJ = diag(wq)*(Vq*J);

while (time<FinalTime)
    
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        timeloc = time + rk4c(INTRK)*dt;
        
        getMesh(meshAcoustic);
        rhsA = AcousRHS2D(Ua,Ue,timeloc);
        getMesh(meshElastic);
        rhsE = ElasRHS2D(Ue,Ua,timeloc);
       
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
        pp(:,(1:Ke)+Ka) = Ue{3}+Ue{4};
        
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

getMesh(meshAcoustic);
err = wqJ.*(Vq*Ua{1} - (u1ax(xq,yq,time) + u2ay(xq,yq,time))).^2;
errA = sqrt(sum(err(:)))

getMesh(meshElastic);
sxx = (2*mu2+lambda) .* u1bx(xq,yq,time) + lambda.*u2by(xq,yq,time);
syy = lambda.*u1bx(xq,yq,time) + (2*mu2+lambda) .* u2by(xq,yq,time);

err = wqJ.*(Vq*(Ue{3}+Ue{4}) - (sxx+syy)).^2;
errE = sqrt(sum(err(:)))


return



function [rhs] = ElasRHS2D(U,Ua,time)

Globals2D;

global Nfld mu lambda Vq Pq tauv taus useWADG
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

opt=3;
if opt==1 % traction BCs
    nSx(mapB) = -2*(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -2*(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
elseif opt==2 % basic ABCs
    nSx(mapB) = -(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
    dU{1}(mapB) = -U{1}(vmapB);
    dU{2}(mapB) = -U{2}(vmapB);
elseif opt==3 % zero velocity
    global v1b v2b
    dU{1}(mapB) = 2*(v1b(x(vmapB),y(vmapB),time)-U{1}(vmapB));
    dU{2}(mapB) = 2*(v2b(x(vmapB),y(vmapB),time)-U{2}(vmapB));
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

flux{1}(mapE) = dSx + tauv(vmapE).*dUn.*nxf;
flux{2}(mapE) = dSy + tauv(vmapE).*dUn.*nyf;
flux{3}(mapE) = dUx.*nxf + taus(vmapE).*(dSx.*nxf);
flux{4}(mapE) = dUy.*nyf + taus(vmapE).*(dSy.*nyf);
flux{5}(mapE) = dUx.*nyf + dUy.*nxf + taus(vmapE).*(dSx.*nyf+dSy.*nxf);

% compute right hand sides of the PDE's
rr{1} =  divSx   +  LIFT*(.5*Fscale.*flux{1});
rr{2} =  divSy   +  LIFT*(.5*Fscale.*flux{2});
rr{3} =  du1dx   +  LIFT*(.5*Fscale.*flux{3});
rr{4} =  du2dy   +  LIFT*(.5*Fscale.*flux{4});
rr{5} =  du12dxy +  LIFT*(.5*Fscale.*flux{5});

rhs{1} = rr{1};
rhs{2} = rr{2};
rhs{3} = (2*mu+lambda).*rr{3} + lambda.*rr{4};
rhs{4} = lambda.*rr{3} + (2*mu+lambda).*rr{4};
rhs{5} = mu .* rr{5};


return;



function rhs = AcousRHS2D(U,Ue,time)

Globals2D;
global cq Vq Pq
global mapA vmapA mapE vmapE

p = U{1}; u = U{2}; v = U{3};

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

ndotdU = nx.*du + ny.*dv;

% % Impose traction
% dp(mapB) = 2*(-p(vmapB));
% ndotdU(mapB) = 0;

% velocity BCs
global v1a v2a
dp(mapB) = 0;
ndotdU(mapB) = 2*((v1a(x(vmapB),y(vmapB),time)-u(vmapB)).*nx(mapB) + ...
    (v2a(x(vmapB),y(vmapB),time)-v(vmapB)).*ny(mapB));

% % ABC
% dp(mapB) = -p(vmapB);
% ndotdU(mapB) = -(nx(mapB).*u(vmapB)+ny(mapB).*v(vmapB));

global tau
fluxp = tau*dp + ndotdU;
fluxu = (tau*ndotdU + dp).*nx;
fluxv = (tau*ndotdU + dp).*ny;

% % compute acoustic-elastic interface flux
nxf = nx(mapA);
nyf = ny(mapA);
Snx = Ue{3}(vmapE).*nxf + Ue{5}(vmapE).*nyf;
Sny = Ue{5}(vmapE).*nxf + Ue{4}(vmapE).*nyf;
dUn = (Ue{1}(vmapE)-u(vmapA)).*nxf + (Ue{2}(vmapE)-v(vmapA)).*nyf;

fluxp(mapA) = tau*((Snx.*nxf + Sny.*nyf) - p(vmapA)) + dUn;
fluxu(mapA) = tau*dUn.*nxf + (Snx - p(vmapA).*nxf);
fluxv(mapA) = tau*dUn.*nyf + (Sny - p(vmapA).*nyf);

pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
rhs{1} =  divU + LIFT*(.5*Fscale.*fluxp);
rhs{2} =  dpdx + LIFT*(.5*Fscale.*fluxu);
rhs{3} =  dpdy + LIFT*(.5*Fscale.*fluxv);

% rhs{1} = Pq*(cq.*(Vq*rhs{1}));
% rhs{1} = c2.*rhs{1};


return

