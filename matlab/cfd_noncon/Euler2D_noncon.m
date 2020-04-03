clear
format long

for N = [2] 
    
    [rqGQ wqGQ] = JacobiGQ(0,0,N);
    [rqGLL wqGLL] = JacobiGL(0,0,N);
    
    sk = 1;
    
    for K1D = [12 24 48]
        fprintf('N = %d, K1D = %d\n',N,K1D)
        
        if ~((N==4) && (K1D == 48))           
            
%             close all
%             a = 0;
%             L2err_GLL{N,sk} = Euler2D(N, K1D, rqGLL, wqGLL, rqGLL, wqGLL,a)

%             close all
%             a = 0;
%             L2err_GQ{N,sk} = Euler2D(N, K1D, rqGQ, wqGQ, rqGQ, wqGQ, a)            
            
            close all
            a = .125;
            L2err_GLL_curved{N,sk} = Euler2D(N, K1D, rqGLL, wqGLL, rqGLL, wqGLL, a)
            
            close all
            a = .125;
            L2err_GQ_curved{N,sk} = Euler2D(N, K1D, rqGQ, wqGQ, rqGQ, wqGQ, a)
            
%             close all
%             a = 0;
%             L2err_mixed{N,sk} = Euler2D(N, K1D, rqGLL, wqGLL, rqGQ, wqGQ, a)
        end
        
        sk = sk + 1;
    end
end

function L2err = Euler2D(Nin, K1D, rq1D, wq1D, rq1D_face, wq1D_face, a)

Globals2D;

if nargin==0
    N = 1;
    K1D = 12;    
    [rq1D wq1D] = JacobiGQ(0,0,N);
    [rq1D_face wq1D_face] = JacobiGQ(0,0,N);    
else
    N = Nin;    
    % rq1D_face = rq1D;    wq1D_face = wq1D;
%     [rq1D_face wq1D_face] = JacobiGQ(0,0,N);    
end
plotMesh = 0;
useNoncon = 1;
if K1D==6
    CFL = .125/2;
else
    CFL = .125;
end
% CFL = 2;
FinalTime = 5;
FinalTime = 1; 
global tau
tau = 1;
% a = .125/2;
% a = .0;
a


% Read in Mesh
if 1
    Lx = 7.5; Ly = 5; ratiox = 1; ratioy = Ly/Lx;
    Kx = 4/3*K1D;
    Ky = K1D;
    [Nv, VX, VY, K, EToV] = QuadMesh2D(Kx,Ky);
else
    Lx = 10; Ly = 5; ratiox = 1; ratioy = Ly/Lx;
    Kx = 2*K1D;
    Ky = K1D;
    [Nv, VX, VY, K, EToV] = QuadMesh2D(Kx,Ky);
end
VX = VX/max(abs(VX));  VY = VY/max(abs(VY)); 
VX = (VX+1)*Lx; VY = VY*Ly;

ids = find(abs(abs(VX)-1) > 1e-8 & abs(abs(VY)-1) > 1e-8);
% VX(ids) = VX(ids) + a*randn(size(ids));
% VY(ids) = VY(ids) + a*randn(size(ids));

StartUpQuad2D;

fmask1   = find( abs(s+1) < NODETOL)';
fmask2   = find( abs(r-1) < NODETOL)';
fmask3   = find( abs(s-1) < NODETOL)'; fmask3 = fmask3(end:-1:1);
fmask4   = find( abs(r+1) < NODETOL)'; fmask4 = fmask4(end:-1:1);
Fmask  = [fmask1;fmask2;fmask3;fmask4]';

BuildPeriodicMaps2D(max(VX)-min(VX),max(VY)-min(VY));

[rp sp] = EquiNodes2D(5); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
% plot(VX,VY,'o')

r1D = JacobiGL(0,0,N);
rp1D = linspace(-1,1,100)';
Vpq1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,rq1D);
Vp1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,r1D);


%%
global M Vq Pq Lq Vff Vf Pfqf VqPq VqLq
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ nrJq nsJq wf
global mapPq


[sq rq] = meshgrid(rq1D);
rq = rq(:); sq = sq(:);
[ws wr] = meshgrid(wq1D);
wq = wr(:).*ws(:);

Vq = Vandermonde2D(N,rq,sq)/V;
M = Vq'*diag(wq)*Vq;

Pq = M\(Vq'*diag(wq)); % J's cancel out

% rq1D_face = rq1D;
% wq1D_face = wq1D;
rq1D = rq1D_face;
wq1D = wq1D_face;

e = ones(size(rq1D));
rf = [rq1D; e; -rq1D; -e];
sf = [-e; rq1D; e; -rq1D];
wf = [wq1D; wq1D; wq1D; wq1D];
Nfq = length(rq1D);

Vf = Vandermonde2D(N,rf,sf)/V;
Lq = M\(Vf'*diag(wf));

nrJ = [0*e; e; 0*e; -e]; % sJ = 2 for all faces,
nsJ = [-e; 0*e; e; 0*e];

Nq = length(rq);
nrJq = repmat(nrJ',Nq,1);
nsJq = repmat(nsJ',Nq,1);

% flux differencing operators
Drq = (Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vf*Pq);
Dsq = (Vq*Ds*Pq - .5*Vq*Lq*diag(nsJ)*Vf*Pq);
VfPq = (Vf*Pq);
VqLq = Vq*Lq;
VqPq = Vq*Pq;

global DNr DNs WN
DNr = [Drq .5*VqLq*diag(nrJ);
    -.5*diag(nrJ)*VfPq .5*diag(nrJ)];
DNs = [Dsq .5*VqLq*diag(nsJ);
    -.5*diag(nsJ)*VfPq .5*diag(nsJ)];
WN = diag([wq;wf]);

QNrskew = .5*(WN*DNr - (WN*DNr)');
QNsskew = .5*(WN*DNs - (WN*DNs)');

DNr = diag(1./[wq;wf])*QNrskew;
DNs = diag(1./[wq;wf])*QNsskew;

DNr = sparse(DNr);
DNs = sparse(DNs);


%% lightweight adaptive mesh data structures

% quadtree
global qt activeK activeK_ids inactiveK_ids Vsplit
qt = cell(K,1);
for e = 1:K
    qt{e} = struct('parent',e,'children',[],'childId',[]);
end

activeK = ones(K,1);

% VDM for split elements and faces
rs = [-1+.5*(1+r); .5*(1+r); .5*(1+r); -1+.5*(1+r)];
ss = [-1+.5*(1+s); -1+.5*(1+s); .5*(1+s); .5*(1+s)];
Vsplit = Vandermonde2D(N,rs,ss)/V;

% split quadrature
% [rq1D wq1D] = JacobiGQ(0,0,N);
rq1Dsplit = [-1+.5*(1+rq1D); .5*(1+rq1D)];
wq1Dsplit = .5*[wq1D; wq1D];

global Tm Efm Emf
Tf = Vandermonde1D(N,rq1D)/Vandermonde1D(N,rq1D_face);
Tm = Vandermonde1D(N,rq1Dsplit)/Vandermonde1D(N,rq1D_face);
Mm = Tm'*diag(wq1Dsplit)*Tm;

% Pm = Mm\(Tm'*diag(wq1Dsplit));
% Pf = Mm\(Tf'*diag(wq1D));

% alternative version: diagonal norm on the faces (for curved)
Mf = Tf'*diag(wq1D)*Tf;  
Pm = Mf\(Tm'*diag(wq1Dsplit));
Pf = Mf\(Tf'*diag(wq1D));

Efm = Tm*Pf; % map from face to mortar nodes
Emf = Tf*Pm; % map from mortar back to face

% plot(Vf*r,Vf*s,'o')
% plot(kron(eye(Nfaces),Efm)*Vf*r,kron(eye(Nfaces),Efm)*Vf*s,'o')
% Emf*Efm
% L2err = [];
% return

% keyboard

%% refine elems

if useNoncon
    Korig = K;
    for ey = 1:Ky
        for ex = 1:Kx
            eid = ex + (ey-1)*Kx + mod(ey,2)-1;
            if mod(ex,2)==0 && eid <= Korig
                hrefine(eid);
            end
        end
    end
else
%    hrefine(1); 
end

activeK_ids = find(activeK);
inactiveK_ids = find(~activeK);


%% build non-con maps

% match faces to faces
ee = {[1 2],[2 3],[3 4],[4 1]}; % new refined elements adjacent to each face
en = {[4 3],[1 4],[2 1],[3 2]}; % neighbors adjacent to each new elem
fn = [3 4 1 2]; % faces of neighbors adjacent to each new elem

% isolate non-conforming faces, find parent elem
split_face_list = zeros(K,Nfaces);
sk = 1;
split_face_counter = 1;
for e = 1:K
    if activeK(e)==1
        for f = 1:Nfaces
            children = qt{EToE(e,f)}.children;
            if isempty(children)
            
            else % if nbr refined, face is non-conforming
                
                ff = f + Nfaces*qt{EToE(e,f)}.children(ee{EToF(e,f)});
                
                split_face_list(e,f) = split_face_counter; 
                split_face_counter = split_face_counter + 1;
            end
        end
        sk = sk + 1;
    end
end

global split_faces
[fid eid] = find(split_face_list'); % transpose for ordering
split_faces = fid + (eid-1)*Nfaces; % list of faces to be split (including non-active elems)
nceid = eid(:);
ncfid = fid(:);
num_conf_faces = Nfaces*K;
num_split_face_list = length(split_faces);
num_faces_total = num_conf_faces+2*num_split_face_list;

EToFnc = zeros(K,Nfaces);
for i = 1:num_split_face_list
    EToFnc(nceid(i),ncfid(i)) = i;
end

xf = reshape(Vf*x,Nfp,Nfaces*K);
yf = reshape(Vf*y,Nfp,Nfaces*K);

% interpolate to mortar points
xfnc = reshape(Tm*xf(:,split_faces),Nfp,num_split_face_list*2);
yfnc = reshape(Tm*yf(:,split_faces),Nfp,num_split_face_list*2);
xf = [xf xfnc];
yf = [yf yfnc];

FToF = zeros(num_faces_total,1); % faces to faces (active faces = nonzero ids)
sk = 1;
for e = 1:K
    if activeK(e) == 1
        for f = 1:Nfaces
            fid = f + (e-1)*Nfaces;
            e_level = get_level(qt,e);
            
            % neighbor data
            enbr = EToE(e,f);
            nbr_children = qt{enbr}.children;
            enbr_level = get_level(qt,enbr);
            
            if ~isempty(nbr_children) % neighbor is refined = split to small face, index into conforming
                
                ncfid = EToFnc(e,f); % get split fid
                enbrs = nbr_children(en{f});
                nbrfid = fn(f) + (enbrs-1)*Nfaces;
                FToF(2*ncfid-1 + num_conf_faces) = nbrfid(1);
                FToF(2*ncfid + num_conf_faces) = nbrfid(2);
                
            elseif e_level==enbr_level+1 % small face to split face: index into split
                
                fnbr = EToFnc(enbr,fn(f)); % neighboring face
                childId = qt{e}.childId;
                if ee{f}(2)==childId % if first part of face f (reversed ordering)
                    FToF(fid) = 2*fnbr-1 + num_conf_faces;
                elseif ee{f}(1)==childId % if second part of face f
                    FToF(fid) = 2*fnbr  + num_conf_faces;
                end
                
            elseif e_level==enbr_level % conforming (elem, nbr same level)
                
                FToF(fid) = EToF(e,f) + (enbr-1)*Nfaces;
                
            else % if all of the above are false, something's wrong!
                
                keyboard
                
            end
        end
        sk = sk + 1;
    end
end

% keyboard

mapMq = reshape(1:Nfp*(num_faces_total),Nfp,num_faces_total);
mapPq = mapMq;

% build node maps using FToF
for f = 1:(num_faces_total)
    fP = FToF(f);
    
    if fP~=0
        
        tol = 1e-10;
        [X1 Y1] = meshgrid(xf(:,f),yf(:,f));
        [X2 Y2] = meshgrid(xf(:,fP),yf(:,fP));
        DX = abs(X1-X2');
        DY = abs(Y1-Y2');
        D = DX + DY;
        [p,~] = find(D<tol);
        
        % NOTE - does not work if K1D too small (1 elem per direction)
        if length(p) == 0
            % assume periodic boundary, find match in x,y
            [px,~] = find(DX<tol);
            [py,~] = find(DY<tol);
            if length(px)==0
                p = py;
            elseif length(py)==0
                p = px;
            else
                keyboard
            end
        end
        
        mapPq(:,f) = mapMq(p,fP);
    end
end

if 0 % check connectivities
    for f = 1:(num_faces_total)
        fP = FToF(f);
        if fP~=0
            for i = 1:Nfp
                id = i + (f-1)*Nfp;
                idM = mapMq(id);
                idP = mapPq(id);
                clf
                plot(xf,yf,'k.','markersize',12)
                hold on
                plot(xf(idM),yf(idM),'o','markersize',16,'linewidth',2)
                plot(xf(idP),yf(idP),'^','markersize',16,'linewidth',2)
                pause%(.1)
            end
        end
    end
end

% return


%% make curvilinear mesh

Lx = max(x(:))-min(x(:));
Ly = max(y(:))-min(y(:));
x = x + a*Lx*cos(pi/Lx.*(x-Lx/2)).*cos(3*pi/Ly.*y);
y = y + a*Ly*sin(4*pi/Lx.*(x-Lx/2)).*cos(pi/Ly.*y);

% % quadratic perturbation
% scale = 15*25;
% dx = (15-x).*(x).*(5-y).*(5+y)/scale;
% dy = (15-x).*(x).*(5-y).*(5+y)/scale;
% x = x + Lx*a*dx;
% y = y + Ly*a*dy;

% inherit nodal positions to make mesh watertight
for e = 1:K
    if length(qt{e}.children)>0 % if e is refined, inherit mapping
        x(:,qt{e}.children) = reshape(Vsplit*x(:,e),Np,4);
        y(:,qt{e}.children) = reshape(Vsplit*y(:,e),Np,4);        
    end
end

xf = reshape(Vf*x,Nfp,Nfaces*K);
yf = reshape(Vf*y,Nfp,Nfaces*K);
xfnc = reshape(Tm*xf(:,split_faces),Nfp,num_split_face_list*2);
yfnc = reshape(Tm*yf(:,split_faces),Nfp,num_split_face_list*2);
xf = [xf xfnc];
yf = [yf yfnc];

xq = Vq*x; yq = Vq*y;
xp = Vp*x; yp = Vp*y;

[rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Vq*Dr,Vq*Ds);
rxJ = rx.*J;    sxJ = sx.*J;
ryJ = ry.*J;    syJ = sy.*J;

[rx,sx,ry,sy,Jf] = GeometricFactors2D(x,y,Vf*Dr,Vf*Ds);
rxJf = rx.*Jf;    sxJf = sx.*Jf;
ryJf = ry.*Jf;    syJf = sy.*Jf;

nxJ = rxJf.*nrJ + sxJf.*nsJ;
nyJ = ryJf.*nrJ + syJf.*nsJ;

sJ = sqrt(nxJ.^2 + nyJ.^2);
nx = nxJ./sJ; ny = nyJ./sJ;

nxJ = reshape(nxJ,Nfp,Nfaces*K);
nyJ = reshape(nyJ,Nfp,Nfaces*K);
sJ = reshape(sJ,Nfp,Nfaces*K);

nxJ = [nxJ reshape(Efm*nxJ(:,split_faces),Nfp,2*num_split_face_list)];
nyJ = [nyJ reshape(Efm*nyJ(:,split_faces),Nfp,2*num_split_face_list)];
sJ = [sJ reshape(Efm*sJ(:,split_faces),Nfp,2*num_split_face_list)];

% interpolate to mortar points
xf = reshape(Vf*x,Nfp,Nfaces*K);
yf = reshape(Vf*y,Nfp,Nfaces*K);
xfnc = reshape(Efm*xf(:,split_faces),Nfp,num_split_face_list*2);
yfnc = reshape(Efm*yf(:,split_faces),Nfp,num_split_face_list*2);
xf = [xf xfnc];
yf = [yf yfnc];

if plotMesh
    
%     plot(xf,yf,'b.','markersize',10)
    hold on
    plot(Vpq1D*xf,Vpq1D*yf,'k-','linewidth',2)
    axis equal
%     for ee = 1:length(activeK_ids)
%         e = activeK_ids(ee);
%         xe = reshape(x(:,e),N+1,N+1);
%         ye = reshape(y(:,e),N+1,N+1);
%         for i = 2:N
%             plot(Vp1D*xe(i,:)',Vp1D*ye(i,:)','k-','linewidth',1)
%         end
%         for j = 2:N
%             plot(Vp1D*xe(:,j),Vp1D*ye(:,j),'k-','linewidth',1)
%         end        
%     end
    axis tight
    axis off
    keyboard
    
%     for ff = 1:length(split_faces)
%         f = split_faces(ff);
%         ncids = Nfaces*K + [2*ff-1 2*ff];
%         
%         plot(Vp1D*xf(:,f),Vp1D*yf(:,f),'b-','linewidth',1) % [xf xfnc],[yf yfnc],'o')        
%         plot(Vp1D*xf(:,ncids),Vp1D*yf(:,ncids),'r--','linewidth',2) % [xf xfnc],[yf yfnc],'o')        
% %         quiver(Tm*xf(:,f),Tm*yf(:,f),Tm*nxJ(:,f),Tm*nyJ(:,f))
% %         quiver(xf(:,f),yf(:,f),nxJ(:,f),nyJ(:,f))
% %         quiver(xf(:,ncids),yf(:,ncids),nxJ(:,ncids),nyJ(:,ncids))
% %         keyboard
% %         pause
%     end
    
    L2err = nan;
    return
end

if 0
    rp1D = linspace(-1,1,100)';
    Vp1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,JacobiGL(0,0,N));
    Vfp = kron(eye(Nfaces),Vp1D);
    xfp = Vfp*x(Fmask(:),:);
    yfp = Vfp*y(Fmask(:),:);
    plot(xfp,yfp,'k.')
    hold on
    plot(x,y,'o')
    L2err = nan;
    return
end

% keyboard

%% fluxes
global gamma
gamma = 1.4;

global V1 V2 V3 V4 U1 U2 U3 U4
rhoe = @(rho,rhou,rhov,E) E - .5*(rhou.^2+rhov.^2)./rho;
pcons = @(rho,rhou,rhov,E) (gamma-1)*rhoe(rho,rhou,rhov,E);
s = @(rho,rhou,rhov,E) log((gamma-1)*rhoe(rho,rhou,rhov,E)./(rho.^gamma));
V1 = @(rho,rhou,rhov,E) (-E + rhoe(rho,rhou,rhov,E).*(gamma + 1 - s(rho,rhou,rhov,E)))./(rhoe(rho,rhou,rhov,E));
V2 = @(rho,rhou,rhov,E) rhou./(rhoe(rho,rhou,rhov,E));
V3 = @(rho,rhou,rhov,E) rhov./(rhoe(rho,rhou,rhov,E));
V4 = @(rho,rhou,rhov,E) (-rho)./(rhoe(rho,rhou,rhov,E));

sV = @(V1,V2,V3,V4) gamma - V1 + (V2.^2+V3.^2)./(2*V4);
rhoeV  = @(V1,V2,V3,V4) ((gamma-1)./((-V4).^gamma)).^(1/(gamma-1)).*exp(-sV(V1,V2,V3,V4)/(gamma-1));
U1 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(-V4);
U2 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(V2);
U3 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(V3);
U4 = @(V1,V2,V3,V4) rhoeV(V1,V2,V3,V4).*(1-(V2.^2+V3.^2)./(2*V4));

global fxS1 fyS1 fxS2 fyS2 fxS3 fyS3 fxS4 fyS4 psix psiy
global pfun beta pavg plogmean vnormavg avg

avg = @(x,y) .5*(x+y);
pfun = @(rho,u,v,E) (gamma-1)*(E-.5*rho.*(u.^2+v.^2));
beta = @(rho,u,v,E) rho./(2*pfun(rho,u,v,E));
pavg     = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER)     avg(rhoL,rhoR)./(2*avg(beta(rhoL,uL,vL,EL),beta(rhoR,uR,vR,ER)));
plogmean = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR)./(2*logmean(beta(rhoL,uL,vL,EL),beta(rhoR,uR,vR,ER)));
vnormavg = @(uL,vL,uR,vR) 2*(avg(uL,uR).^2 + avg(vL,vR).^2) - (avg(uL.^2,uR.^2) + avg(vL.^2,vR.^2));

% entropy potentials
psix = @(rho,u,v,E) (gamma-1)*rho.*u;
psiy = @(rho,u,v,E) (gamma-1)*rho.*v;


%% problem params setup

[rhoq uq vq pq] = vortexSolution(xq,yq,0);
rho  = VqPq*rhoq;
rhou = VqPq*(rhoq.*uq);
rhov = VqPq*(rhoq.*vq);
E    = VqPq*(pq/(gamma-1) + .5*rhoq.*(uq.^2+vq.^2));


%%

global wJq
wJq = diag(wq)*(J);

% Runge-Kutta residual storage
res1 = zeros(Nq,K);
res2 = zeros(Nq,K);
res3 = zeros(Nq,K);
res4 = zeros(Nq,K);

% compute time step size
CN = (N+1)*(N+2); % constant in trace ineq
% CNh = max(CN*max(sJ(:))./min(Jf(:)));
CNh = CN*K1D/(max(VY)-min(VY)); % estimate CFL
dt = CFL*2/CNh;

Nsteps = ceil(FinalTime/dt)
dt = FinalTime/Nsteps;

entropy = zeros(Nsteps,1);
figure(1)
for i = 1:Nsteps
    for INTRK = 1:5
        
        % project to entropy variables for curvilinear
        q1 = Pq*V1(rho,rhou,rhov,E);
        q2 = Pq*V2(rho,rhou,rhov,E);
        q3 = Pq*V3(rho,rhou,rhov,E);
        q4 = Pq*V4(rho,rhou,rhov,E);
        
        % evaluate at quad/surface points
        q1q = Vq*q1;
        q2q = Vq*q2;
        q3q = Vq*q3;
        q4q = Vq*q4;
        
        q1M = Vf*q1;
        q2M = Vf*q2;
        q3M = Vf*q3;
        q4M = Vf*q4;
        
        rhoq  = U1(q1q,q2q,q3q,q4q);
        uq    = U2(q1q,q2q,q3q,q4q)./rhoq;
        vq    = U3(q1q,q2q,q3q,q4q)./rhoq;
        Eq    = U4(q1q,q2q,q3q,q4q);
        
        [rhs1 rhs2 rhs3 rhs4]  = RHS2Dsimple(rhoq,uq,vq,Eq,q1M,q2M,q3M,q4M);

        if (INTRK==5)
            tmp = wJq.*(q1q.*rhs1 + q2q.*rhs2 + q3q.*rhs3 + q4q.*rhs4);
            rhstest(i) = sum(sum(tmp(:,activeK_ids)));
        end
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        res3 = rk4a(INTRK)*res3 + dt*rhs3;
        res4 = rk4a(INTRK)*res4 + dt*rhs4;
        
        rho  = rho  + rk4b(INTRK)*res1;
        rhou = rhou + rk4b(INTRK)*res2;
        rhov = rhov + rk4b(INTRK)*res3;
        E    = E    + rk4b(INTRK)*res4;
        
    end;
    
    Sq = -rho.*s(rho,rhou,rhov,E);
    entropy(i) = sum(sum(wJq(:,activeK_ids).*Sq(:,activeK_ids)));
    
    if  (mod(i,10)==0 || i==1 || i==Nsteps)
        clf
        pp = rho(:,activeK_ids);
        vv = real(Vp*Pq*pp);
        color_line3(xp(:,activeK_ids),yp(:,activeK_ids),vv,vv,'.');
        axis equal
        axis tight
        title(sprintf('time = %f, step %d out of %d, N = %d, K1D = %d, rhstest = %g',dt*i,i,Nsteps,N,K1D,rhstest(i)))
        drawnow
    end
end


[rq2 sq2 wq2] = Cubature2D(Nq+4);
Vq2 = Vandermonde2D(N,rq2,sq2)/V;
xq2 = Vq2*x; yq2 = Vq2*y;
wJq2 = diag(wq2)*(Vq2*Pq*J);
[rhoex uex vex pex] = vortexSolution(xq2,yq2,FinalTime);

rhouex = rhoex.*uex;
rhovex = rhoex.*vex;
% p = (gamma-1)*(E-.5*rho*(u^2+v^2));
Eex = pex/(gamma-1) + .5*rhoex.*(uex.^2+vex.^2);

rhoq = Vq2*Pq*rho;
rhouq = Vq2*Pq*rhou;
rhovq = Vq2*Pq*rhov;
Eq = Vq2*Pq*E;
% err = wJq2.*((rhoq-rhoex).^2 + (rhouq-rhouex).^2 + (rhovq-rhovex).^2 + (Eq-Eex).^2);
% err = err(:,activeK_ids);
% L2err = sqrt(sum(err(:)));

rhoerr = max(max(abs(rhoq(:,activeK_ids)-rhoex(:,activeK_ids))));
rhouerr = max(max(abs(rhouq(:,activeK_ids)-rhouex(:,activeK_ids))));
rhoverr = max(max(abs(rhovq(:,activeK_ids)-rhovex(:,activeK_ids))));
Eerr = max(max(abs(rhovq(:,activeK_ids)-rhovex(:,activeK_ids))));
L2err = max([rhoerr,rhouerr,rhoverr,Eerr]); % hack - L2 = Linf error here

dS = abs(entropy-entropy(1));
dS = entropy + max(abs(entropy));
figure(2)
semilogy(dt*(1:Nsteps),dS,'o--');hold on
semilogy(dt*(1:Nsteps),abs(rhstest),'x--');hold on
legend('dS','rhstest')

end

function [rhs1 rhs2 rhs3 rhs4] = RHS2Dsimple(rhoq,uq,vq,Eq,v1f,v2f,v3f,v4f)

Globals2D;

global M Vq Pq Lq Lqf Vff Vf VqPq VqLq
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ nrJq nsJq
global mapPq
global fxS1 fyS1 fxS2 fyS2 fxS3 fyS3 fxS4 fyS4
global pfun beta pavg plogmean vnormavg avg
global Drq Dsq VfPq

global DNr DNs

%% ------------------
% projection stuff

global Efm Emf
global split_faces
global V1 V2 V3 V4 U1 U2 U3 U4
global gamma tau

v1f = reshape(v1f,Nfp,Nfaces*K);
v2f = reshape(v2f,Nfp,Nfaces*K);
v3f = reshape(v3f,Nfp,Nfaces*K);
v4f = reshape(v4f,Nfp,Nfaces*K);

v1fnc = [v1f reshape(Efm*v1f(:,split_faces),Nfp,2*length(split_faces))];
v2fnc = [v2f reshape(Efm*v2f(:,split_faces),Nfp,2*length(split_faces))];
v3fnc = [v3f reshape(Efm*v3f(:,split_faces),Nfp,2*length(split_faces))];
v4fnc = [v4f reshape(Efm*v4f(:,split_faces),Nfp,2*length(split_faces))];

rhoM = U1(v1fnc,v2fnc,v3fnc,v4fnc);
uM   = U2(v1fnc,v2fnc,v3fnc,v4fnc)./rhoM;
vM   = U3(v1fnc,v2fnc,v3fnc,v4fnc)./rhoM;
EM   = U4(v1fnc,v2fnc,v3fnc,v4fnc);
betaM = beta(rhoM,uM,vM,EM);

rhoP = rhoM(mapPq);
uP = uM(mapPq);
vP = vM(mapPq);
EP = EM(mapPq);
betaP = betaM(mapPq);

[FxS1 FxS2 FxS3 FxS4 FyS1 FyS2 FyS3 FyS4] = euler_flux(rhoM,rhoP,uM,uP,vM,vP,betaM,betaP,gamma);
fSf1 = nxJ.*FxS1 + nyJ.*FyS1;
fSf2 = nxJ.*FxS2 + nyJ.*FyS2;
fSf3 = nxJ.*FxS3 + nyJ.*FyS3;
fSf4 = nxJ.*FxS4 + nyJ.*FyS4;

% Lax-Friedrichs flux
unorm2 = (uM.^2+vM.^2);
pM = (gamma-1)*(EM - .5*rhoM.*unorm2);
cvel = sqrt(gamma*pM./rhoM);
lam = sqrt(unorm2)+cvel;
LFc = max(lam(mapPq),lam).*sJ;
fSf1 = fSf1 - .25*tau*LFc.*(rhoP-rhoM);
fSf2 = fSf2 - .25*tau*LFc.*(rhoP.*uP-rhoM.*uM);
fSf3 = fSf3 - .25*tau*LFc.*(rhoP.*vP-rhoM.*vM);
fSf4 = fSf4 - .25*tau*LFc.*(EP-EM);

for ff = 1:length(split_faces)
    f = split_faces(ff);
    ncids = Nfaces*K + [2*ff-1 2*ff];
    
    % flux vs mortar points
    rhom = reshape(rhoM(:,ncids),2*Nfp,1);
    um = reshape(uM(:,ncids),2*Nfp,1);
    vm = reshape(vM(:,ncids),2*Nfp,1);
    betam = reshape(betaM(:,ncids),2*Nfp,1);    
    
    [rhox, rhoy] = meshgrid(rhom,rhoM(:,f));
    [ux, uy] = meshgrid(um,uM(:,f));
    [vx, vy] = meshgrid(vm,vM(:,f));
    [betax, betay] = meshgrid(betam,betaM(:,f));
    [rhoxlog, rhoylog] = meshgrid(log(rhom),log(rhoM(:,f)));
    [betaxlog, betaylog] = meshgrid(log(betam),log(betaM(:,f)));
    
    [nxJ1, nxJ2] = meshgrid(reshape(nxJ(:,ncids),2*Nfp,1),nxJ(:,f));
    [nyJ1, nyJ2] = meshgrid(reshape(nyJ(:,ncids),2*Nfp,1),nyJ(:,f));
    nxJa = .5*(nxJ1+nxJ2);
    nyJa = .5*(nyJ1+nyJ2);
    
    % flux correction terms: E
    [FxS1 FxS2 FxS3 FxS4 FyS1 FyS2 FyS3 FyS4] = euler_flux(rhox,rhoy,ux,uy,vx,vy,betax,betay,gamma,...
        rhoxlog,rhoylog,betaxlog,betaylog);
    
    fSc1 = nxJa.*FxS1 + nyJa.*FyS1;
    fSc2 = nxJa.*FxS2 + nyJa.*FyS2;
    fSc3 = nxJa.*FxS3 + nyJa.*FyS3;
    fSc4 = nxJa.*FxS4 + nyJa.*FyS4;
    
    fc1 = sum(Emf.*fSc1,2) - Emf*sum(Efm.*fSc1',2);
    fc2 = sum(Emf.*fSc2,2) - Emf*sum(Efm.*fSc2',2);
    fc3 = sum(Emf.*fSc3,2) - Emf*sum(Efm.*fSc3',2);
    fc4 = sum(Emf.*fSc4,2) - Emf*sum(Efm.*fSc4',2);
    
    fSf1(:,f) = Emf*reshape(fSf1(:,ncids),2*Nfp,1) + fc1;
    fSf2(:,f) = Emf*reshape(fSf2(:,ncids),2*Nfp,1) + fc2;
    fSf3(:,f) = Emf*reshape(fSf3(:,ncids),2*Nfp,1) + fc3;
    fSf4(:,f) = Emf*reshape(fSf4(:,ncids),2*Nfp,1) + fc4;
end

% get rid of extra face info
fSf1 = reshape(fSf1(:,1:Nfaces*K),Nfp*Nfaces,K);
fSf2 = reshape(fSf2(:,1:Nfaces*K),Nfp*Nfaces,K);
fSf3 = reshape(fSf3(:,1:Nfaces*K),Nfp*Nfaces,K);
fSf4 = reshape(fSf4(:,1:Nfaces*K),Nfp*Nfaces,K);

rhoM = reshape(rhoM(:,1:Nfaces*K),Nfp*Nfaces,K);
uM = reshape(uM(:,1:Nfaces*K),Nfp*Nfaces,K);
vM = reshape(vM(:,1:Nfaces*K),Nfp*Nfaces,K);
EM = reshape(EM(:,1:Nfaces*K),Nfp*Nfaces,K);
betaM = reshape(betaM(:,1:Nfaces*K),Nfp*Nfaces,K);

%%
global activeK_ids inactiveK_ids

% combine vol/flux vars
rho = [rhoq; rhoM];
u = [uq; uM];
v = [vq; vM];
betaN = [beta(rhoq,uq,vq,Eq);betaM];

divF1 = zeros(size(DNr,1),K);
divF2 = zeros(size(DNr,1),K);
divF3 = zeros(size(DNr,1),K);
divF4 = zeros(size(DNr,1),K);
for ee = 1:length(activeK_ids)
    
    e = activeK_ids(ee);
    
    [rhox, rhoy] = meshgrid(rho(:,e));
    [ux, uy] = meshgrid(u(:,e));
    [vx, vy] = meshgrid(v(:,e));    
    [betax, betay] = meshgrid(betaN(:,e));
    [rhoxlog, rhoylog] = meshgrid(log(rho(:,e)));
    [betaxlog, betaylog] = meshgrid(log(betaN(:,e)));
    [FxS1 FxS2 FxS3 FxS4 FyS1 FyS2 FyS3 FyS4] = euler_flux(rhox,rhoy,ux,uy,vx,vy,betax,betay,gamma,...
        rhoxlog,rhoylog,betaxlog,betaylog);
    
    [rxJ1, rxJ2] = meshgrid([rxJ(:,e);rxJf(:,e)]);
    [sxJ1, sxJ2] = meshgrid([sxJ(:,e);sxJf(:,e)]);
    [ryJ1, ryJ2] = meshgrid([ryJ(:,e);ryJf(:,e)]);
    [syJ1, syJ2] = meshgrid([syJ(:,e);syJf(:,e)]);
    rxJK = avg(rxJ1,rxJ2);  sxJK = avg(sxJ1,sxJ2);
    ryJK = avg(ryJ1,ryJ2);  syJK = avg(syJ1,syJ2);
    
    Dx = DNr.*rxJK + DNs.*sxJK;
    Dy = DNr.*ryJK + DNs.*syJK;
    
    divF1(:,e) = sum(Dx.*FxS1,2) + sum(Dy.*FyS1,2);
    divF2(:,e) = sum(Dx.*FxS2,2) + sum(Dy.*FyS2,2);
    divF3(:,e) = sum(Dx.*FxS3,2) + sum(Dy.*FyS3,2);
    divF4(:,e) = sum(Dx.*FxS4,2) + sum(Dy.*FyS4,2);
end

PN = 2*[VqPq VqLq];
rhs1 = PN*divF1 + VqLq*(fSf1);
rhs2 = PN*divF2 + VqLq*(fSf2);
rhs3 = PN*divF3 + VqLq*(fSf3);
rhs4 = PN*divF4 + VqLq*(fSf4);

% assume diag mass
rhs1 = -(rhs1./J);
rhs2 = -(rhs2./J);
rhs3 = -(rhs3./J);
rhs4 = -(rhs4./J);

end


function [FxS1 FxS2 FxS3 FxS4 FyS1 FyS2 FyS3 FyS4] = ...
    euler_flux(rhox,rhoy,ux,uy,vx,vy,betax,betay,gamma,rhoxlog,rhoylog,betaxlog,betaylog)

if nargin == 9
    rholog = logmean(rhox,rhoy);
    betalog = logmean(betax,betay);
else
    rholog = logmean(rhox,rhoy,rhoxlog,rhoylog);
    betalog = logmean(betax,betay,betaxlog,betaylog);
end

% arithmetic avgs
rhoavg = .5*(rhox+rhoy);
uavg = .5*(ux+uy);
vavg = .5*(vx+vy);

vnavg = ux.*uy + vx.*vy; % 2*(uavg.^2 + vavg.^2) - .5*((ux.^2+uy.^2) + (vx.^2+vy.^2));
pa = rhoavg./(betax+betay);
f4aux = rholog./(2*(gamma-1)*betalog) + pa + .5*rholog.*vnavg;

FxS1 = rholog.*uavg;       
FxS2 = FxS1.*uavg + pa;    
FxS3 = FxS1.*vavg;               
FxS4 = f4aux.*uavg;

FyS1 = rholog.*vavg; 
FyS2 = FxS3; %FyS1.*uavg;
FyS3 = FyS1.*vavg + pa;
FyS4 = f4aux.*vavg;

end

function [rho u v p] = vortexSolution(x,y,t)

global gamma
x0 = 5;
y0 = 0;
beta = 5;
r2 = (x-x0-t).^2 + (y-y0).^2;

u = 1 - beta*exp(1-r2).*(y-y0)/(2*pi);
v = beta*exp(1-r2).*(x-x0-t)/(2*pi);
rho = 1 - (1/(8*gamma*pi^2))*(gamma-1)/2*(beta*exp(1-r2)).^2;
rho = rho.^(1/(gamma-1));
p = rho.^gamma;

% rho = (2 + 0*sin(pi*(x - t)));
% u = ones(size(x));
% v = zeros(size(x));
% p = ones(size(x));

if 0
    % pulse condition
    x0 = 5;
    rho = 2 + 0*(abs(x-x0) < 5);
    u = rho;
    v = 0*rho;
    p = rho.^gamma;
    
end

end

%% other routines: refining 1-irregular meshes

%         nbr 3
%      2----7----3
%      | e4 | e3 |
% nbr4 8----9----6  nbr 2
%      | e1 | e2 |
%      1----5----4
%         nbr 1
function hrefine(e)

Globals2D
global qt activeK Vsplit Tm

if activeK(e)==0
    fprintf('Cannot refine element %d: not active\n',e);
    return
end

% TODO: enforce 1-irregularity

% add to quadtree
qt = [qt;cell(4,1)];
qt{e}.children = K+(1:4);
qt{K+1} = struct('parent',e,'children',[],'childId',1);
qt{K+2} = struct('parent',e,'children',[],'childId',2);
qt{K+3} = struct('parent',e,'children',[],'childId',3);
qt{K+4} = struct('parent',e,'children',[],'childId',4);

% update neighbor data
EToE = [EToE;zeros(4,Nfaces)]; % expand neighbor array
EToF = [EToF;zeros(4,Nfaces)]; % expand neighbor array
newnbr = K + (1:4);

% local neighbors
EToE(K+1,2) = newnbr(2); EToE(K+1,3) = newnbr(4);
EToE(K+2,4) = newnbr(1); EToE(K+2,3) = newnbr(3);
EToE(K+3,1) = newnbr(2); EToE(K+3,4) = newnbr(4);
EToE(K+4,1) = newnbr(1); EToE(K+4,2) = newnbr(3);
%         nbr 3
%      2----7----3
%      | e4 | e3 |
% nbr4 8----9----6  nbr 2
%      | e1 | e2 |
%      1----5----4
%         nbr 1
EToF(K+1,2) = 4; EToF(K+1,3) = 1;
EToF(K+2,4) = 2; EToF(K+2,3) = 1;
EToF(K+3,1) = 3; EToF(K+3,4) = 2;
EToF(K+4,1) = 3; EToF(K+4,2) = 4;

% inherit external neighbors: check nbrs of e
ee = {[1 2],[2 3],[3 4],[4 1]}; % new refined elements adjacent to each face
en = {[4 3],[1 4],[2 1],[3 2]}; % neighbors adjacent to each new elem
fn = [3 4 1 2]; % faces of neighbors adjacent to each new elem
for f = 1:4
    for i = 1:2 % neighbors per face
        new_elem = K+ee{f}(i);
        if isempty(qt{EToE(e,f)}.children)
            % if neighbor is not refined, connect to parent
            parent = qt{EToE(e,f)}.parent;
            EToE(new_elem,f) = parent;
            EToF(new_elem,f) = EToF(parent,f);
        else % if neighbor refined, pick up neighbors from children
            nbr_elem = qt{EToE(e,f)}.children(en{f}(i));
            EToE(new_elem,f) = nbr_elem;
            EToE(nbr_elem,fn(f)) = new_elem;
            
            EToF(new_elem,f) = EToF(e,f); % inherit face connectivity
        end
    end
end

% could also build by sweeping quadtree but this is easier
activeK(e) = 0;
activeK = [activeK; ones(4,1)];

xs = reshape(Vsplit*x(:,e),Np,4);
ys = reshape(Vsplit*y(:,e),Np,4);
x = [x xs];
y = [y ys];

% update number of elements
K = K + 4; % count refined elements here

end

function level = get_level(qt,e)

next_elem = e;
level = 0;
at_top = 0;
while ~at_top
    parent = qt{next_elem}.parent;
    if parent == next_elem
        at_top = 1;
    else
        next_elem = parent;
        level = level + 1;
    end
end
end
