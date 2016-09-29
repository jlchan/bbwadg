function EulerVortex2D

clear -global *
clear

Globals2D;

N = 7;
nref = 1;
useJprojection = 1;

Nq = 2*N+1;

filename = 'Grid/Other/circA01.neu';
filename = 'Grid/CNS2D/couette.neu';
filename = 'Grid/CFD/pvortexA025.neu';
% filename = 'Grid/Other/block2.neu';
[Nv, VX, VY, K, EToV, BCType] = MeshReaderGambitBC2D(filename);
VX = VX/max(abs(VX));  VY = VY/max(abs(VY));
VX = VX*5 + 5; VY = VY*5;

% K1D = 2;
% [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

StartUp2D;
BuildBCMaps2D;
for ref = 1:nref
    Refine2D(ones(size(EToV)));
    StartUp2D;
    BuildBCMaps2D;
end

% BuildPeriodicMaps2D(10,10);
% PlotMesh2D;axis on;return

if 0
    Bfaces = zeros(size(EToV));
    for e = 1:K
        Bfaces(e,EToE(e,:)==e) = 1;
        rad = sqrt(VX(EToV(e,:)).^2 + VY(EToV(e,:)).^2);
        if any(rad < .5)
            Bfaces(e,EToE(e,:)==e) = 2;
        end
    end
    [k1,f1] = find(Bfaces==1);
    [k2,f2] = find(Bfaces==2);
    MakeCylinder2D([k1,f1], 1, 0, 0);
    MakeCylinder2D([k2,f2], 1/4, 0, 0); % double cylinder    
end
cinfo = BuildCurvedOPS2D_opt(Nq,useJprojection);

% PlotMesh2D; axis on;return

% jesse custom for curved driver
global Vq wq Vrq Vsq 
global Vfq wfq VfqFace xf yf wq
global rxq sxq ryq syq Jq Jfq nxq nyq sJq 
global Pq Prq Psq Pfq 

[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

[rq sq wq] = Cubature2D(Nq);
% [rq sq wq] = QNodes2D(N); [rq sq] = xytors(rq,sq);
[rq1D, wq1D] = JacobiGQ(0, 0, Nq);
rfq = [rq1D; -rq1D; -ones(size(rq1D))];
sfq = [-ones(size(rq1D)); rq1D; -rq1D];
wfq = repmat(wq1D,Nfaces,1);

Vq = Vandermonde2D(N,rq,sq)/V;
[Vrq Vsq] = GradVandermonde2D(N,rq,sq);
Vrq = Vrq/V; Vsq = Vsq/V;
Pq = (V*V')*(Vq'*diag(wq));
Prq = (V*V')*(Vrq'*diag(wq));
Psq = (V*V')*(Vsq'*diag(wq));

global Vcolloc Drc Dsc
[rc sc wc] = QNodes2D(N); [rc sc] = xytors(rc,sc);
[rc sc] = Nodes2D(N); [rc sc] = xytors(rc,sc);
% rc = [-0.8968532316898645,0.3982993518653447,0.3978382631466154,-0.4980940100997226,-0.8998412094465894,0.7920309320167708,-0.9001058622488827,-0.05007732972494768,0.02289750084338348,-0.8979020577966027,-0.05161942712715207,-0.4975235348941438,-0.5110877092127721,-0.896172961623547,-0.5119017247820049]';
% sc = [-0.0530686701240277,-0.4984580098292808,-0.9003147500848969,0.3982002632311986,-0.4984583976019817,-0.8941299128326118,0.3982009522673672,-0.05306872288989844,-0.5109958715674856,-0.8941298666552433,-0.8967623483705202,-0.9003148085904928,0.02217564219171834,0.7923450404031387,-0.5109965442976414]';
% rc = [-0.577481890524152,-0.6979784209072049,-0.9261005378863045,0.1568879701435784,-0.9045760697551653,-0.9027097501052324,0.5920624018508198,-0.5047786755466918,-0.6238190954657445,-0.05811159673672549,0.8360908598546688,-0.2713883345507647,-0.6368305132454305,-0.9230064512583414,0.07720573184150727,0.5467234608560792,-0.2956492619400888,-0.9150781568290514,0.2359040658118374,-0.923864588768103,-0.3098475385273428]';
% sc = [0.5004883552652438,-0.8940839938189025,-0.15110519726495,-0.6521092849831685,0.1759645278878392,-0.9333811241802712,-0.8940839939774403,-0.6521092997468331,0.2476381922174764,-0.3050578769882315,-0.9333811176159436,0.1759644108911185,-0.3050578905992031,0.5004883823376606,-0.1511051853246766,-0.6316453075555499,-0.408701458747774,-0.6316453169671435,-0.9260565364813174,0.8477291791681128,-0.9260565325714376]';

Vc = Vandermonde2D(N,rc,sc);
Vcolloc = Vc/V; % interp from nodal to colloc pts
[Vrc Vsc] = GradVandermonde2D(N,rc,sc);
Drc = Vrc/Vc; Dsc = Vsc/Vc; % collocation derivative pts

% interp nodal to quadrature
Vfq = Vandermonde2D(N,rfq,sfq)/V;

Nfq = size(Vfq,1)/Nfaces;
xf = Vfq*x; yf = Vfq*y;
[mapM,mapP] = BuildNodeMaps2D(xf,yf,EToE,EToF);
mapM = reshape(mapM,Nfq*Nfaces,K);
mapP = reshape(mapP,Nfq*Nfaces,K);

if 1
    % make periodic box BCs
    mapBf = find(mapM==mapP);
    tol = 1e-8;
    dX = max(VX)-min(VX);
    dY = max(VY)-min(VY);
    mapPp = mapP;
    for i = 1:length(mapP(:))
        id1 = mapM(i);
        id2 = mapP(i);
        if (id1==id2) % boundary node
            dx = abs(abs(xf(id1)-xf(mapBf))-dX) + abs(yf(id1)-yf(mapBf));
            dy = abs(xf(id1)-xf(mapBf)) + abs(abs(yf(id1)-yf(mapBf))-dY);
            if min(dx(:)) < tol
                mapPp(i) = mapBf(dx(:)<tol);
            elseif min(dy(:))<tol
                mapPp(i) = mapBf(dy(:)<tol);
            end
        end
    end
    mapP = mapPp; % set periodic maps
end

V1D = Vandermonde1D(N,JacobiGL(0,0,N));
VfqFace = Vandermonde1D(N,rq1D)/V1D;
VfqFace = blkdiag(VfqFace,VfqFace,VfqFace); % repeat for 3 faces
Pfq = (V*V')*(Vfq'*diag(wfq));

Nc = length(wq); Nfc = length(cinfo(1).gnx(:));
rxq = zeros(Nc,K); sxq = zeros(Nc,K); 
ryq = zeros(Nc,K); syq = zeros(Nc,K); 
Jq = zeros(Nc,K);
nxq = zeros(Nfc,K); nyq = zeros(Nfc,K);
sJq = zeros(Nfc,K); Jfq = zeros(Nfc,K); 
for e = 1:K
    [rxk,sxk,ryk,syk,Jk] = GeometricFactors2D(x(:,e),y(:,e),Vrq,Vsq);
    rxq(:,e) = rxk; sxq(:,e) = sxk; 
    ryq(:,e) = ryk; syq(:,e) = syk; 
    Jq(:,e) = Jk; 
    
    nxq(:,e) = cinfo(e).gnx(:); nyq(:,e) = cinfo(e).gny(:);
    sJq(:,e) = cinfo(e).gsJ(:); Jfq(:,e) = cinfo(e).gJ(:);    
end

global tau;
tau = 1;

global gamma
gamma = 1.4;

%%

if 0
    for t = 0:.1:5
        [rho u v p] = vortexSolution(xp,yp,t);
        m1 = rho.*u;
        m2 = rho.*v;
        E = p/(gamma-1) + .5*rho.*(u.^2 + v.^2);
        
        vv = v;
        clf
        color_line3(xp,yp,vv,vv,'.')
        view(3)
        drawnow
    end
    return
end

%%

% Set initial conditions
d = 5; x0 = 4.9; y0 = 2.45;
rho = exp(-d*sqrt((x-x0).^2 + (y-y0).^2));
m1 = zeros(Np, K); m2 = zeros(Np, K); E = zeros(Np,K);

[rho, u, v, p] = vortexSolution(x,y,0);
m1 = rho.*u;
m2 = rho.*v;
E = p/(gamma-1) + .5*rho.*(u.^2 + v.^2);

FinalTime = 4;

% setup
resrho = zeros(Np,K);resm1 = zeros(Np,K); resm2 = zeros(Np,K); resE = zeros(Np,K);

rLGL = JacobiGL(0,0,N); rmin = abs(rLGL(2)-rLGL(1));
dt = .1 * rmin / max(Fscale(:));

[re se] = EquiNodes2D(N); [re se] = xytors(re,se);
Ve = Vandermonde2D(N,re,se)/V;
xe = Ve*x; ye = Ve*y;
VB = bern_basis_tri(N,r,s);

% outer time step loop
time = 0; tstep = 0;
while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
        
    for INTRK = 1:5
%         [rhsm1, rhsm2, rhsrho] = acousticsRHS2D_JC(m1,m2,rho); rhsE = 0*rhsm1;
%         [rhsrho, rhsm1, rhsm2, rhsE] = EulerRHS2D_strong(rho,m1,m2,E);
%         [rhsrho, rhsm1, rhsm2, rhsE] = EulerRHS2D_colloc(rho,m1,m2,E);
        [rhsrho, rhsm1, rhsm2, rhsE] = EulerRHS2D(rho,m1,m2,E);       
        
        % initiate and increment Runge-Kutta residuals
        resrho = rk4a(INTRK)*resrho + dt*rhsrho;
        resm1 = rk4a(INTRK)*resm1 + dt*rhsm1;
        resm2 = rk4a(INTRK)*resm2 + dt*rhsm2;
        resE = rk4a(INTRK)*resE + dt*rhsE;
        
        % update fields
        rho = rho + rk4b(INTRK)*resrho;
        m1  = m1 + rk4b(INTRK)*resm1;
        m2  = m2 + rk4b(INTRK)*resm2;
        E   = E + rk4b(INTRK)*resE;        
    end
    
    if 1 && mod(tstep,10)==0
        clf
        pp = Vp*rho;
%         pp = Vp*m2;
        color_line3(xp,yp,pp,pp,'.');
        hold on
%         rhoB = VB\rho;
%         plot3(xe,ye,rhoB,'o');
%         [TV TVx TVy] = TV2D(rhoB);
%         TV = repmat(TV,size(Vp,1),1);
%         color_line3(xp,yp,TV,TV,'.')
%         view(90,0)
%         axis([-1 1 -1 1 -.5 .5])
        title(sprintf('time = %f',time));
        drawnow
    end
    
    % Increment time
    time = time+dt; 
    tstep = tstep+1;
    if mod(tstep,10) ==0
        disp(sprintf('on timestep %d out of %d\n',tstep,round(FinalTime/dt)))
    end
end




function [rhsm1, rhsm2, rhsrho] = acousticsRHS2D_JC(u,v,p)

% function [rhsm1, rhsm2, rhsrho] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

global Vq Vrq Vsq VfqFace
global rxq sxq ryq syq nxq nyq  
global Pq Pfq Jq sJq

Globals2D;

% Define field differences at faces
pf = VfqFace*p(Fmask(:),:);
uf = VfqFace*u(Fmask(:),:);
vf = VfqFace*v(Fmask(:),:);
Nfq = size(VfqFace,1)/Nfaces;
dp = zeros(Nfq*Nfaces,K); dp(:) = pf(mapP)-pf(mapM);
du = zeros(Nfq*Nfaces,K); du(:) = uf(mapP)-uf(mapM);
dv = zeros(Nfq*Nfaces,K); dv(:) = vf(mapP)-vf(mapM);

% % Impose reflective boundary conditions (p+ = -p-)
% du(mapB) = 0; dv(mapB) = 0; dp(mapB) = -2*p(vmapB);

% % can interp numerical fluxes *after* since mult by sJ is higher order
% dp = VfqFace*dp;
% du = VfqFace*du;
% dv = VfqFace*dv;

% evaluate upwind fluxes
ndotdU = nxq.*du + nyq.*dv;
global tau
fluxp =  tau*dp - ndotdU;
fluxu =  (tau*ndotdU - dp).*nxq;
fluxv =  (tau*ndotdU - dp).*nyq;

% local derivatives of fields
pr = Vrq*p; ps = Vsq*p;
ur = Vrq*u; us = Vsq*u;
vr = Vrq*v; vs = Vsq*v;

dpdx = pr.*rxq + ps.*sxq;
dpdy = pr.*ryq + ps.*syq;
dudx = ur.*rxq + us.*sxq;
dvdy = vr.*ryq + vs.*syq;
divU = dudx + dvdy; 

% % compute right hand sides of the PDE's
% % Petrov-Galerkin DG
% rhsrho =  Pq*(-divU) + Pfq*(fluxp.*sJq./Jfq/2.0);
% rhsm1 =  Pq*(-dpdx) + Pfq*(fluxu.*sJq./Jfq/2.0);
% rhsm2 =  Pq*(-dpdy) + Pfq*(fluxv.*sJq./Jfq/2.0);

% compute right hand sides of the PDE's
rhsrho =  Pq*(-divU.*Jq) + Pfq*(fluxp.*sJq/2.0);
rhsm1 =  Pq*(-dpdx.*Jq) + Pfq*(fluxu.*sJq/2.0);
rhsm2 =  Pq*(-dpdy.*Jq) + Pfq*(fluxv.*sJq/2.0);

% apply inverse mass matrix
rhsrho = Pq*((Vq*rhsrho)./Jq);
rhsm1 = Pq*((Vq*rhsm1)./Jq);
rhsm2 = Pq*((Vq*rhsm2)./Jq);

return;


function [F G rho u v p] = EulerFluxes(rho,m1,m2,E)

global gamma

u = m1./rho; 
v = m2./rho; 
p = (gamma-1)*(E - .5*rho.*(u.^2 + v.^2));

F{1} = rho.*u;
F{2} = rho.*u.^2 + p;
F{3} = rho.*u.*v;
F{4} = u.*(E+p);

G{1} = rho.*v;
G{2} = rho.*u.*v;
G{3} = rho.*v.^2 + p;
G{4} = v.*(E+p);



function [rhsrho, rhsm1, rhsm2, rhsE] = EulerRHS2D_colloc(rho,m1,m2,E)

global gamma
global Vq Vrq Vsq VfqFace Vfq xf yf wq
global rxq sxq ryq syq nxq nyq 
global Pq Prq Psq Pfq Jq sJq
global Vcolloc Drc Dsc

Globals2D;

rhoq = Vcolloc*rho; 
m1q  = Vcolloc*m1; 
m2q  = Vcolloc*m2; 
Eq   = Vcolloc*E;
[Fq, Gq] = EulerFluxes(rhoq,m1q,m2q,Eq);
for fld = 1:4
    Fq{fld} = Vcolloc\Fq{fld};
    Gq{fld} = Vcolloc\Gq{fld};
end

Qf{1} = Vfq * rho; 
Qf{2} = Vfq * m1; 
Qf{3} = Vfq * m2; 
Qf{4} = Vfq * E;  
[F, G, rhof, uf, vf, pf] = EulerFluxes(Qf{1},Qf{2},Qf{3},Qf{4});

Nfq = size(rhof,1)/Nfaces;
for fld = 1:4            
    
    FM = F{fld}(mapM);   FP = F{fld}(mapP);     
    GM = G{fld}(mapM);   GP = G{fld}(mapP);
    QM = Qf{fld}(mapM);  QP = Qf{fld}(mapP);
    
    qjump = QM - QP;
    Favg = .5*(FP-FM);
    Gavg = .5*(GP-GM);
    
    % Lax Friedrichs flux
    lambdaM = sqrt(uf(mapM).^2 + vf(mapM).^2) + sqrt(gamma*abs(pf(mapM)./rhof(mapM)));
    lambdaP = sqrt(uf(mapP).^2 + vf(mapP).^2) + sqrt(gamma*abs(pf(mapP)./rhof(mapP)));
    lambda = zeros(Nfq*Nfaces,K);
    lambda(:) = max(lambdaM(:),lambdaP(:));
    for f = 1:Nfaces
        ids = 1:Nfq;
        lambda(ids + (f-1)*Nfq,:) = repmat(max(lambda(ids + (f-1)*Nfq,:),[],1),Nfq,1);
    end
    dF = Favg.*nxq + Gavg.*nyq + .5*lambda.*qjump; % roe flux
    
    Fqr = Dr * Fq{fld};    Fqs = Ds * Fq{fld};
    Gqr = Dr * Gq{fld};    Gqs = Ds * Gq{fld};
    dFdx = rx.*Fqr + sx.*Fqs;
    dGdy = ry.*Gqr + sy.*Gqs;
    divF = (dFdx + dGdy).*J;
    
    rhs{fld} = -divF - Pfq*(dF.*sJq);
end

% curvilinear correction
rhsrho = rhs{1}./J;
rhsm1 = rhs{2}./J;
rhsm2 = rhs{3}./J;
rhsE = rhs{4}./J;

return;

function [rhsrho, rhsm1, rhsm2, rhsE] = EulerRHS2D_strong(rho,m1,m2,E)

global gamma
global Vq Vrq Vsq VfqFace Vfq xf yf wq
global rxq sxq ryq syq nxq nyq 
global Pq Prq Psq Pfq Jq sJq

Globals2D;

rhoq = Vq*rho; 
m1q  = Vq*m1; 
m2q  = Vq*m2; 
Eq   = Vq*E;
[Fq, Gq] = EulerFluxes(rhoq,m1q,m2q,Eq);

for fld = 1:4        
    Fq{fld} = Pq*Fq{fld};
    Gq{fld} = Pq*Gq{fld};
end

Qf{1} = Vfq * rho; 
Qf{2} = Vfq * m1; 
Qf{3} = Vfq * m2; 
Qf{4} = Vfq * E;  
[F, G, rhof, uf, vf, pf] = EulerFluxes(Qf{1},Qf{2},Qf{3},Qf{4});

Nfq = size(rhof,1)/Nfaces;
for fld = 1:4            
    
    FM = F{fld}(mapM);   FP = F{fld}(mapP);     
    GM = G{fld}(mapM);   GP = G{fld}(mapP);
    QM = Qf{fld}(mapM);  QP = Qf{fld}(mapP);
    
    qjump = QM - QP;
    Favg = .5*(FP-FM);
    Gavg = .5*(GP-GM);
    
    % Lax Friedrichs flux
    lambdaM = sqrt(uf(mapM).^2 + vf(mapM).^2) + sqrt(gamma*abs(pf(mapM)./rhof(mapM)));
    lambdaP = sqrt(uf(mapP).^2 + vf(mapP).^2) + sqrt(gamma*abs(pf(mapP)./rhof(mapP)));
    lambda = zeros(Nfq*Nfaces,K);
    lambda(:) = max(lambdaM(:),lambdaP(:));
    for f = 1:Nfaces
        ids = 1:Nfq;
        lambda(ids + (f-1)*Nfq,:) = repmat(max(lambda(ids + (f-1)*Nfq,:),[],1),Nfq,1);
    end
    dF = Favg.*nxq + Gavg.*nyq + .5*lambda.*qjump; % roe flux
    
    Fqr = Dr * Fq{fld};    Fqs = Ds * Fq{fld};
    Gqr = Dr * Gq{fld};    Gqs = Ds * Gq{fld};
    dFdx = rx.*Fqr + sx.*Fqs;
    dGdy = ry.*Gqr + sy.*Gqs;
    divF = (dFdx + dGdy).*J;
    
    rhs{fld} = -divF - Pfq*(dF.*sJq);
end

% curvilinear correction
rhsrho = rhs{1}./J;
rhsm1 = rhs{2}./J;
rhsm2 = rhs{3}./J;
rhsE = rhs{4}./J;

return;


function [rhsrho, rhsm1, rhsm2, rhsE] = EulerRHS2D(rho,m1,m2,E)

global gamma
global Vq Vrq Vsq VfqFace Vfq xf yf
global rxq sxq ryq syq nxq nyq 
global Pq Prq Psq Pfq Jq sJq

Globals2D;

rhoq = Vq*rho; 
m1q  = Vq*m1; 
m2q  = Vq*m2; 
Eq   = Vq*E;
[Fq, Gq] = EulerFluxes(rhoq,m1q,m2q,Eq);

Qf{1} = VfqFace * rho(Fmask(:),:);
Qf{2} = VfqFace * m1(Fmask(:),:);
Qf{3} = VfqFace * m2(Fmask(:),:);
Qf{4} = VfqFace * E(Fmask(:),:);
[F, G, rhof, uf, vf, pf] = EulerFluxes(Qf{1},Qf{2},Qf{3},Qf{4});

% inflow_ids = find(abs(xf(mapP))<1e-4 | abs(abs(yf(mapP))-5)<1e-4);
% [rhoBC uBC vBC pBC] = vortexSolution(xf(mapP),yf(mapP),time);
% m1BC = rhoBC.*uBC;
% m2BC = rhoBC.*vBC;
% EBC = pBC/(gamma-1) + .5*rhoBC.*(uBC.^2 + vBC.^2);
% [FBC GBC] = EulerFluxes(rhoBC,m1BC,m2BC,EBC);
% QBC{1} = rhoBC;
% QBC{2} = m1BC;
% QBC{3} = m2BC;
% QBC{4} = EBC;

Nfq = size(rhof,1)/Nfaces;
for fld = 1:4            
    
    FM = F{fld}(mapM);   FP = F{fld}(mapP);     
    GM = G{fld}(mapM);   GP = G{fld}(mapP);
    QM = Qf{fld}(mapM);  QP = Qf{fld}(mapP);
                
%     FP(inflow_ids) = FBC{fld}(inflow_ids);
%     GP(inflow_ids) = GBC{fld}(inflow_ids);
%     QP(inflow_ids) = QBC{fld}(inflow_ids);
    
    qjump = QM - QP;
    Favg = FM + FP;
    Gavg = GM + GP;
    
    % Lax Friedrichs flux    
    lambdaM = sqrt(uf(mapM).^2 + vf(mapM).^2) + sqrt(gamma*abs(pf(mapM)./rhof(mapM)));
    lambdaP = sqrt(uf(mapP).^2 + vf(mapP).^2) + sqrt(gamma*abs(pf(mapP)./rhof(mapP)));
    lambda = zeros(Nfq*Nfaces,K);
    lambda(:) = max(lambdaM(:),lambdaP(:));
    for f = 1:Nfaces
        ids = 1:Nfq;
        lambda(ids + (f-1)*Nfq,:) = repmat(max(lambda(ids + (f-1)*Nfq,:),[],1),Nfq,1);
    end
    dF = Favg.*nxq + Gavg.*nyq + lambda.*qjump; % roe flux
    
    FGr = (rxq.*Fq{fld} + ryq.*Gq{fld});
    FGs = (sxq.*Fq{fld} + syq.*Gq{fld}); 
    divF = Prq*(FGr.*Jq) + Psq*(FGs.*Jq);
    
    rhs{fld} = divF - Pfq*(dF.*sJq/2);
end

% curvilinear correction
rhsrho = Pq*((Vq*rhs{1})./Jq);
rhsm1 = Pq*((Vq*rhs{2})./Jq);
rhsm2 = Pq*((Vq*rhs{3})./Jq);
rhsE = Pq*((Vq*rhs{4})./Jq);

return;

function [rho u v p] = vortexSolution(x,y,t)

global gamma
x0 = 4; 
y0 = 0;
beta = 5;
r = sqrt((x-x0-t).^2 + (y-y0).^2);

u = 1 - beta*exp(1-r.^2).*(y-y0)/(2*pi);
v = beta*exp(1-r.^2).*(x-x0-t)/(2*pi);
rho = 1 - (gamma-1)/(16*gamma*pi^2)*beta^2*exp(2*(1-r.^2));
rho = rho.^(1/(gamma-1));
p = rho.^gamma;


% copied in b/c matlab is funny
function MakeCylinder2D(faces, ra,xo,yo)

% Function: MakeCylinder2D(faces, ra, xo, yo)
% Purpose:  Use Gordon-Hall blending with an isoparametric map to modify a list
%           of faces so they conform to a cylinder of radius r centered at (xo,yo)
Globals2D;

NCurveFaces = size(faces,1);
vflag = zeros(size(VX));
for n=1:NCurveFaces
    
    % move vertices of faces to be curved onto circle
    k = faces(n,1); f = faces(n,2);
    v1 = EToV(k, f); v2 = EToV(k, mod(f,Nfaces)+1);
    
    % compute polar angles of start and end face vertices relative to circle center
    theta1 = atan2(VY(v1)-yo,VX(v1)-xo);
    theta2 = atan2(VY(v2)-yo,VX(v2)-xo);
    
    % move vertices onto circle
    newx1 = xo + ra*cos(theta1); newy1 = yo + ra*sin(theta1);
    newx2 = xo + ra*cos(theta2); newy2 = yo + ra*sin(theta2);
    
    % update mesh vertex locations
    VX(v1) = newx1; VX(v2) = newx2; VY(v1) = newy1; VY(v2) = newy2;
    
    % store modified vertex numbers
    vflag(v1) = 1;  vflag(v2) = 1;
end

% map modified vertex flag to each element
vflag = vflag(EToV);

% locate elements with at least one modified vertex
ks = find(sum(vflag,2)>0);

% build coordinates of all the corrected nodes
va = EToV(ks,1)'; vb = EToV(ks,2)'; vc = EToV(ks,3)';
x(:,ks) = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));
y(:,ks) = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));

for n=1:NCurveFaces  % deform specified faces
    k = faces(n,1); f = faces(n,2);
    
    % find vertex locations for this face and tangential coordinate
    if(f==1) v1 = EToV(k,1); v2 = EToV(k,2); vr = r; end
    if(f==2) v1 = EToV(k,2); v2 = EToV(k,3); vr = s; end
    if(f==3) v1 = EToV(k,1); v2 = EToV(k,3); vr = s; end
    fr = vr(Fmask(:,f));
    x1 = VX(v1); y1 = VY(v1); x2 = VX(v2); y2 = VY(v2);
    
    % move vertices at end points of this face to the cylinder
    theta1 = atan2(y1-yo, x1-xo); theta2 = atan2(y2-yo, x2-xo);
    
    % check to make sure they are in the same quadrant
    if ((theta2 > 0) & (theta1 < 0)), theta1 = theta1 + 2*pi; end;
    if ((theta1 > 0) & (theta2 < 0)), theta2 = theta2 + 2*pi; end;
    
    % distribute N+1 nodes by arc-length along edge
    theta = 0.5*theta1*(1-fr) + 0.5*theta2*(1+fr);

    % evaluate deformation of coordinates
    fdx = xo + ra*cos(theta)-x(Fmask(:,f),k);
    fdy = yo + ra*sin(theta)-y(Fmask(:,f),k);

    % build 1D Vandermonde matrix for face nodes and volume nodes
    Vface = Vandermonde1D(N, fr);  Vvol  = Vandermonde1D(N, vr);
    % compute unblended volume deformations
    vdx = Vvol*(Vface\fdx); vdy = Vvol*(Vface\fdy);

    % blend deformation and increment node coordinates
    ids = find(abs(1-vr)>1e-7); % warp and blend
    if(f==1) blend = -(r(ids)+s(ids))./(1-vr(ids)); end;
    if(f==2) blend =      +(r(ids)+1)./(1-vr(ids)); end;
    if(f==3) blend = -(r(ids)+s(ids))./(1-vr(ids)); end;
    
    x(ids,k) = x(ids,k)+blend.*vdx(ids);
    y(ids,k) = y(ids,k)+blend.*vdy(ids);

end

% repair other coordinate dependent information
Fx = x(Fmask(:), :); Fy = y(Fmask(:), :);
[rx,sx,ry,sy,J] = GeometricFactors2D(x, y,Dr,Ds);
[nx, ny, sJ] = Normals2D(); Fscale = sJ./(J(Fmask,:));
return

function [cinfo] = BuildCurvedOPS2D_opt(intN,useJprojection)

% function [cinfo] = BuildCurvedOPS2D(intN)
% Purpose: build curved info for curvilinear elements
Globals2D;

% 1. Create cubature information
% 1.1 Extract cubature nodes and weights
[cR,cS,cW,Ncub] = Cubature2D(intN);

% 1.1. Build interpolation matrix (nodes->cubature nodes)
cV = InterpMatrix2D(cR, cS);

% 1.2 Evaluate derivatives of Lagrange interpolants at cubature nodes
[cDr,cDs] = Dmatrices2D(N, cR, cS, V);

% 2. Create surface quadrature information

% 2.1 Compute Gauss nodes a nd weights for 1D integrals
[gz, gw] = JacobiGQ(0, 0, intN);

% 2.2 Build Gauss nodes running counter-clockwise on element faces
gR = [gz, -gz, -ones(size(gz))];
gS = [-ones(size(gz)), gz, -gz];

% 2.3 For each face
for f1=1:Nfaces
    % 2.3.1 build nodes->Gauss quadrature interpolation and differentiation matrices
    gV(:,:,f1) = InterpMatrix2D(gR(:,f1), gS(:,f1));
    [gDr(:,:,f1),gDs(:,:,f1)] = Dmatrices2D(N, gR(:,f1), gS(:,f1), V);
end

% 3. For each curved element, evaluate custom operator matrices
[k,f] = find(BCType);
% curved = sort(unique(k));
curved = 1:K; % force all elems = curved
Ncurved = length(curved);

% 3.1 Store custom information in array of Matlab structs
cinfo = [];
for c=1:Ncurved
    % find next curved element and the coordinates of its nodes
    k1 = curved(c); x1 = x(:,k1); y1 = y(:,k1); cinfo(c).elmt = k1;
    
    % compute geometric factors
    [crx,csx,cry,csy,cJ] = GeometricFactors2D(x1,y1,cDr,cDs);
    
    % build mass matrix
    if useJprojection==1 % project using quadrature
        if (c==1)
            disp('using WADG')
        end
        
        M = cV'*diag(cW)*cV;
        MinvJ = cV'*diag(cW./cJ)*cV;
        cMM = M*(MinvJ\M);

    else % true mass matrix
        if (c==1)
            disp('using true mass')
        end
        cMM = cV'*diag(cJ.*cW)*cV;
        
    end
    cinfo(c).MM = cMM;
    
    % build physical derivative matrices
    cinfo(c).Dx = cMM\(cV'*diag(cW.*cJ)*(diag(crx)*cDr + diag(csx)*cDs));
    cinfo(c).Dy = cMM\(cV'*diag(cW.*cJ)*(diag(cry)*cDr + diag(csy)*cDs));
    
    % build individual lift matrices at each face
    for f1=1:Nfaces
        k2 = EToE(k1,f1); f2 = EToF(k1,f1);
        
        % compute geometric factors
        [grx,gsx,gry,gsy,gJ] = GeometricFactors2D(x1,y1,gDr(:,:,f1),gDs(:,:,f1));
        
        % compute normals and surface Jacobian at Gauss points on face f1
        if(f1==1) gnx = -gsx;     gny = -gsy;     end;
        if(f1==2) gnx =  grx+gsx; gny =  gry+gsy; end;
        if(f1==3) gnx = -grx;     gny = -gry;     end;
        
        gsJ = sqrt(gnx.*gnx+gny.*gny);       
        gnx = gnx./gsJ;  gny = gny./gsJ;  gsJ = gsJ.*gJ;
        
        cinfo(c).gsJ(:,f1) = gsJ;
        cinfo(c).gJ(:,f1) = gJ;
        
        % store normals and coordinates at Gauss nodes
        cinfo(c).gnx(:,f1) = gnx;  cinfo(c).gx(:,f1)  = gV(:,:,f1)*x1;
        cinfo(c).gny(:,f1) = gny;  cinfo(c).gy(:,f1)  = gV(:,:,f1)*y1;
        
        % store Vandermondes for '-' and '+' traces
        cinfo(c).gVM(:,:,f1) = gV(:,:,f1);
        cinfo(c).gVP(:,:,f1) = gV(end:-1:1,:,f2);
        
        % compute and store matrix to lift Gauss node data
        cinfo(c).glift(:,:,f1) = cMM\(gV(:,:,f1)'*diag(gw.*gsJ));
    end
end

return