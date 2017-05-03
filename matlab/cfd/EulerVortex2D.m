function EulerVortex2D

clear -global *
clear

Globals2D;

N = 4;

FinalTime = 15;

Nq = 2*N+1;

K1D = 4;
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
VX = VX/max(abs(VX));  VY = VY/max(abs(VY));
VX = VX*5 + 5; VY = VY*5;

StartUp2D;

% jesse custom for curved driver
global Vq wq Vrq Vsq 
global Vfq wfq VfqFace xf yf 
global rxq sxq ryq syq Jq Jfq nxq nyq sJq 
global Pq Prq Psq Pfq 

[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

[rq sq wq] = Cubature2D(Nq);
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
% [rc sc] = Nodes2D(N); [rc sc] = xytors(rc,sc);
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

Nc = length(wq); Nfc = size(VfqFace,1);
rxq = zeros(Nc,K); sxq = zeros(Nc,K); 
ryq = zeros(Nc,K); syq = zeros(Nc,K); 
Jq = zeros(Nc,K);
nxq = zeros(Nfc,K); nyq = zeros(Nfc,K);
sJq = zeros(Nfc,K); Jfq = zeros(Nfc,K); 

[rxq,sxq,ryq,syq,Jq] = GeometricFactors2D(x,y,Vrq,Vsq);
    
nxq = VfqFace*nx; nyq = VfqFace*ny;
sJq = VfqFace*sJ; Jfq = VfqFace*J(Fmask(:),:);

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
        
        vv = E;
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

% setup
resrho = zeros(Np,K);resm1 = zeros(Np,K); resm2 = zeros(Np,K); resE = zeros(Np,K);

rLGL = JacobiGL(0,0,N); rmin = abs(rLGL(2)-rLGL(1));
dt = .1 * rmin / max(Fscale(:));

[re se] = EquiNodes2D(N); [re se] = xytors(re,se);
Ve = Vandermonde2D(N,re,se)/V;
xe = Ve*x; ye = Ve*y;
VB = bern_basis_tri(N,r,s);

M = inv(V*V');

f = zeros(Np,1);
sk = 1;
for i = 0:N
    for j = 0:N-i
        f(sk) = 1;
        if (i==N || j==N)
            f(sk) = .1;
        end
        sk = sk + 1;
    end
end
F = V * diag(f) *inv(V);

% keyboard

% outer time step loop
time = 0; tstep = 1;
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
        
        % filter
        rho = F*rho;
        m1 = F*m1;
        m2 = F*m2;
        E = F*E;
                
    end
    
    enorm = J.*(M*(rho.^2 + m1.^2 + m2.^2 + E.^2));    
    energy(tstep) = sqrt(sum(sum(enorm)));
    
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
        colorbar
        drawnow
    end
    
    % Increment time
    time = time+dt; 
    tstep = tstep+1;
    if mod(tstep,10) ==0
        disp(sprintf('on timestep %d out of %d\n',tstep,round(FinalTime/dt)))
    end
end

figure
plot(energy)




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

% curvilinear correction using collocation
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
%rho = 1 - (gamma-1)/(16*gamma*pi^2)*beta^2*exp(2*(1-r.^2));
rho = 1 - (1/(8*gamma*pi^2))*(gamma-1)/2*(beta*exp(1-r.^2)).^2;
rho = rho.^(1/(gamma-1));
p = rho.^gamma;



