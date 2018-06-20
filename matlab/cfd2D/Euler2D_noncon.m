clear -globals
clear

% useQuads = 1;mypath

Globals2D;
N = 1;
K1D = 2;
useSkew = 1;
CFL = .25;
FinalTime = 1;
global tau
tau = 1;

% Read in Mesh
Lx = 10; Ly = 5; ratiox = 1; ratioy = Ly/Lx;
%[Nv, VX, VY, K, EToV] = unif_tri_mesh(round(ratiox*K1D),round(K1D*ratioy));
[Nv, VX, VY, K, EToV] = QuadMesh2D(round(ratiox*K1D),round(K1D*ratioy));
VX = VX/max(abs(VX));  VY = VY/max(abs(VY));
VX = (VX+1)*Lx; VY = VY*Ly;

% [Nv, VX, VY, K, EToV] = QuadMesh2D(K1D,K1D);

a = 1/K1D;
wadgProjEntropyVars = abs(a)>1e-8;

ids = find(abs(abs(VX)-1) > 1e-8 & abs(abs(VY)-1) > 1e-8);
% VX(ids) = VX(ids) + a*randn(size(ids));
% VY(ids) = VY(ids) + a*randn(size(ids));

StartUpQuad2D;

% PlotMesh2D
% for e = 1:K
%     vids = EToV(e,:);
%     text(mean(VX(vids)),mean(VY(vids)),num2str(e))
% end
% return

fmask1   = find( abs(s+1) < NODETOL)'; 
fmask2   = find( abs(r-1) < NODETOL)';
fmask3   = find( abs(s-1) < NODETOL)'; fmask3 = fmask3(end:-1:1);
fmask4   = find( abs(r+1) < NODETOL)'; fmask4 = fmask4(end:-1:1);
Fmask  = [fmask1;fmask2;fmask3;fmask4]';

BuildPeriodicMaps2D(max(VX)-min(VX),max(VY)-min(VY));

[rp sp] = EquiNodes2D(15); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

%
global M Vq Pq Lq Vfqf Vfq Pfqf VqPq VqLq
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ nrJq nsJq wfq
global mapPq

% define vol quad pts
[rq1D wq1D] = JacobiGQ(0,0,N);
[sq rq] = meshgrid(rq1D);
rq = rq(:); sq = sq(:);
[ws wr] = meshgrid(wq1D); 
wq = wr(:).*ws(:);

Vq = Vandermonde2D(N,rq,sq)/V;
M = Vq'*diag(wq)*Vq;

Pq = M\(Vq'*diag(wq)); % J's cancel out
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

rxJ = rx.*J; sxJ = sx.*J;
ryJ = ry.*J; syJ = sy.*J;
rxJ = Vq*rxJ; sxJ = Vq*sxJ;
ryJ = Vq*ryJ; syJ = Vq*syJ;
J = Vq*J;

[rq1D wq1D] = JacobiGQ(0,0,N);
% [rq1D wq1D] = JacobiGQ(0,0,N); % exact for 2N-1
% [rq1D2 wq1D2] = JacobiGQ(0,0,N);
% [rq1D wq1D] = JacobiGQ(0,0,N-1); % exact for 2N-1
% [rq1D wq1D] = JacobiGQ(0,0,N-1); % exact for 2N-1
% VV = Vandermonde1D(N,rq1D); 
% rq1D = [-1+.5*(1+rq1D); .5*(1+rq1D)]; % split quadrature
% wq1D = .5*[wq1D;wq1D];

e = ones(size(rq1D));
rfq = [rq1D; e; -rq1D; -e]; sfq = [-e; rq1D; e; -rq1D]; wfq = [wq1D; wq1D; wq1D; wq1D];
Vfqf = kron(eye(Nfaces),Vandermonde1D(N,rq1D)/Vandermonde1D(N,JacobiGL(0,0,N)));
% rfq = [rq1D2; e; -rq1D2; -e]; sfq = [-e; rq1D; e; -rq1D]; wfq = [wq1D2; wq1D; wq1D2; wq1D];
% Vfqf = blkdiag(Vandermonde1D(N,rq1D2)/Vandermonde1D(N,JacobiGL(0,0,N)),Vandermonde1D(N,rq1D)/Vandermonde1D(N,JacobiGL(0,0,N)),...
%     Vandermonde1D(N,rq1D2)/Vandermonde1D(N,JacobiGL(0,0,N)),Vandermonde1D(N,rq1D)/Vandermonde1D(N,JacobiGL(0,0,N)));
Nfq = length(rq1D);
% plot(rfq,sfq,'o');return


Vfq = Vandermonde2D(N,rfq,sfq)/V;

Mf = Vfq'*diag(wfq)*Vfq;
Lq = M\(Vfq'*diag(wfq));


nx = Vfqf*nx;
ny = Vfqf*ny;
sJ = Vfqf*sJ;
nxJ = (nx.*sJ);
nyJ = (ny.*sJ);
Fscale = Vfqf*Fscale;

nrJ = [0*e; e; 0*e; -e]; % sJ = 2 for all faces, 
nsJ = [-e; 0*e; e; 0*e];
% quiver(rfq,sfq,nrJ,nsJ);return

Nq = length(rq);
nrJq = repmat(nrJ',Nq,1);
nsJq = repmat(nsJ',Nq,1);

% flux differencing operators

Drq = (Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vfq*Pq);
Dsq = (Vq*Ds*Pq - .5*Vq*Lq*diag(nsJ)*Vfq*Pq);
VfPq = (Vfq*Pq);
VqLq = Vq*Lq;
VqPq = Vq*Pq;

global DNr DNs WN 
DNr = [Vq*Dr*Pq-.5*Vq*Lq*diag(nrJ)*Vfq*Pq .5*VqLq*diag(nrJ);
    -.5*diag(nrJ)*VfPq .5*diag(nrJ)];
DNs = [(Vq*Ds*Pq-.5*Vq*Lq*diag(nsJ)*Vfq*Pq) .5*VqLq*diag(nsJ);
    -.5*diag(nsJ)*VfPq .5*diag(nsJ)];
WN = diag([wq;wfq]);

% weak derivative operators
Drqw = (Vq*(M\(Dr'*Vq'*diag(wq))) - .5*Vq*Lq*diag(nrJ)*Vfq*Pq);
Dsqw = (Vq*(M\(Ds'*Vq'*diag(wq))) - .5*Vq*Lq*diag(nsJ)*Vfq*Pq);

% DNrw = [Drqw -.5*VqLq*diag(nrJ);
%     .5*diag(nrJ)*VfPq .5*diag(nrJ)];
% DNsw = [Dsqw -.5*VqLq*diag(nsJ);
%     .5*diag(nsJ)*VfPq .5*diag(nsJ)];

QNrskew = .5*(WN*DNr - (WN*DNr)');
QNsskew = .5*(WN*DNs - (WN*DNs)');

if useSkew
    DNr = diag(1./[wq;wfq])*QNrskew;
    DNs = diag(1./[wq;wfq])*QNsskew;    
end

%% lightweight adaptive mesh data structures

global qt activeK Vsplit Vfsplit

% quadtree
qt = cell(K,1);
for e = 1:K
    qt{e} = struct('parent',e,'children',[],'childId',[]);
end

activeK = ones(K,1);

% VDM for split elements and faces
rs = [-1+.5*(1+r); .5*(1+r); .5*(1+r); -1+.5*(1+r)];
ss = [-1+.5*(1+s); -1+.5*(1+s); .5*(1+s); .5*(1+s)];
Vsplit = Vandermonde2D(N,rs,ss)/V;

rq1Dsplit = [-1+.5*(1+rq1D); .5*(1+rq1D)];
wq1Dsplit = .5*[wq1D; wq1D];
Vfsplit = Vandermonde1D(N,rq1Dsplit)/Vandermonde1D(N,rq1D);
Pqsplit = (Vfsplit'*diag(wq1Dsplit)*Vfsplit)\(Vfsplit'*diag(wq1Dsplit));

%% refine elems

%hrefine(round(K/2-K1D/2));
% hrefine(1);
% hrefine(6);
% hrefine(7);
% hrefine(10);
% hrefine(11);
% hrefine(10);


%% build non-con maps

% isolate non-conforming faces, find parent elem           
nonconFaces = zeros(K,Nfaces);
sk = 1;
noncon_face_counter = 1;
for e = 1:K
    if activeK(e)==1
        for f = 1:Nfaces
            children = qt{EToE(e,f)}.children;
            if isempty(children) 
            else % if nbr refined, face is non-conforming                
                nonconFaces(e,f) = noncon_face_counter; %noncon_face_counter;
                noncon_face_counter = noncon_face_counter + 1;
            end
        end        
        sk = sk + 1;
    end
end

global ncfaces
[fid eid] = find(nonconFaces'); % transpose for ordering
ncfaces = fid + (eid-1)*Nfaces; % list of faces to be split (including non-active elems)
nceid = eid(:); 
ncfid = fid(:); 
num_conf_faces = Nfaces*K;
num_nonconf_faces = length(ncfaces);
num_faces_total = num_conf_faces+2*num_nonconf_faces;

EToFnc = zeros(K,Nfaces);
for i = 1:num_nonconf_faces
    EToFnc(nceid(i),ncfid(i)) = i;
end

% xf = reshape(x(Fmask(:),:),Nfp,Nfaces*K);
% yf = reshape(y(Fmask(:),:),Nfp,Nfaces*K);
xf = reshape(Vfq*x,Nfp,Nfaces*K);
yf = reshape(Vfq*y,Nfp,Nfaces*K);
xfnc = reshape(Vfsplit*xf(:,ncfaces),Nfp,num_nonconf_faces*2);
yfnc = reshape(Vfsplit*yf(:,ncfaces),Nfp,num_nonconf_faces*2);

eK = find(activeK);
% text(mean(x(:,eK)),mean(y(:,eK)),num2str(eK))
% hold on
% plot([xf xfnc],[yf yfnc],'o')
% xf(:,ncfaces) = Pqsplit*reshape(xfnc,2*Nfp,num_nonconf_faces);
% yf(:,ncfaces) = Pqsplit*reshape(yfnc,2*Nfp,num_nonconf_faces);
% plot(xf,yf,'^')
% return

% match faces to faces
ee = {[1 2],[2 3],[3 4],[4 1]}; % new refined elements adjacent to each face
en = {[4 3],[1 4],[2 1],[3 2]}; % neighbors adjacent to each new elem
fn = [3 4 1 2]; % faces of neighbors adjacent to each new elem

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
                              
                
                ncfid = EToFnc(e,f); % get noncon fid
                enbrs = nbr_children(en{f});
                nbrfid = fn(f) + (enbrs-1)*Nfaces;
                FToF(2*ncfid-1 + num_conf_faces) = nbrfid(1);
                FToF(2*ncfid + num_conf_faces) = nbrfid(2);

            elseif e_level==enbr_level+1 % small face to split face: index into noncon
                                              
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

xf = [xf xfnc];
yf = [yf yfnc];

mapMq = reshape(1:Nfp*(num_faces_total),Nfp,num_faces_total);
mapPq = mapMq;

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
                
        % NOTE - does not work if K1D is too small!!
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

% for f = 1:(num_faces_total)
%     fP = FToF(f);   
%     if fP~=0 
%         for i = 1:Nfp
%             id = i + (f-1)*Nfp;
%             idM = mapMq(id);
%             idP = mapPq(id);
%             clf
%             plot(xf,yf,'k.','markersize',12)
%             hold on
%             plot(xf(idM),yf(idM),'o','markersize',16,'linewidth',2)
%             plot(xf(idP),yf(idP),'^','markersize',16,'linewidth',2)
%             pause(.1)
%         end
%     end
% end
% return


%% make curvilinear mesh 

x0 = Lx; y0 = 0;
% x0 = 0; y0 = 0; Lx = 1; Ly = 1;
x = x + Lx*a*cos(1/2*pi*(x-x0)/Lx).*cos(3/2*pi*(y-y0)/Ly);
y = y + Ly*a*sin(pi*(x-x0)/Lx).*cos(1/2*pi*(y-y0)/Ly);

xq = Vq*x; yq = Vq*y;
xp = Vp*x; yp = Vp*y;
xf = Vfq*x;    yf = Vfq*y;

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

[rxk,sxk,ryk,syk,J] = GeometricFactorsQuad2D(x,y,Vq*Dr,Vq*Ds);
rxJ = rxk.*J;
ryJ = ryk.*J;
sxJ = sxk.*J;
syJ = syk.*J;

[rxk,sxk,ryk,syk,Jk] = GeometricFactorsQuad2D(x,y,Vfq*Dr,Vfq*Ds);
rxJf = rxk.*Jk;    sxJf = sxk.*Jk;
ryJf = ryk.*Jk;    syJf = syk.*Jk;
Jf = Jk;

for e = 1:K
%     [rxk,sxk,ryk,syk,Jk] = GeometricFactorsQuad2D(x(:,e),y(:,e),Vq*Dr,Vq*Ds);
%     rxJ(:,e) = rxk.*Jk;    sxJ(:,e) = sxk.*Jk;
%     ryJ(:,e) = ryk.*Jk;    syJ(:,e) = syk.*Jk;
%     J(:,e) = Jk;
    
%     [rxk,sxk,ryk,syk,Jk] = GeometricFactorsQuad2D(x(:,e),y(:,e),Vfq*Dr,Vfq*Ds);
%     rxJf(:,e) = rxk.*Jk;    sxJf(:,e) = sxk.*Jk;
%     ryJf(:,e) = ryk.*Jk;    syJf(:,e) = syk.*Jk;
%     Jf(:,e) = Jk;
end

nxJ = rxJf.*nrJ + sxJf.*nsJ;
nyJ = ryJf.*nrJ + syJf.*nsJ;

sJ = sqrt(nxJ.^2 + nyJ.^2);
nx = nxJ./sJ; ny = nyJ./sJ;
return
% nx = nxJ./Jf;
% ny = nyJ./Jf;
% sJ = sqrt(nx.^2 + ny.^2);
% nx = nx./sJ; ny = ny./sJ;
% sJ = sJ.*Jf;

% e = 1;
% DNx = .5*(diag([rxJ(:,e); rxJf(:,e)])*DNr + DNr*diag([rxJ(:,e); rxJf(:,e)]) ...
%      + diag([sxJ(:,e); sxJf(:,e)])*DNs + DNs*diag([sxJ(:,e); sxJf(:,e)]));
% BNr = [zeros(Nq) zeros(Nq,Nfq*Nfaces);
%     zeros(Nfq*Nfaces,Nq) diag(wfq.*nrJ)];
% BNs = [zeros(Nq) zeros(Nq,Nfq*Nfaces);
%     zeros(Nfq*Nfaces,Nq) diag(wfq.*nsJ)];
% 
% BNx = diag([rxJ(:,e); rxJf(:,e)])*BNr + diag([sxJ(:,e); sxJf(:,e)])*BNs;
% QNx = WN*DNx;
% 
% u = randn(size(DNx,2),1);
% e = ones(size(DNx,2),1);
% e'*QNx*u
% e'*(BNx-QNx')*u

% keyboard

%% fluxes
global gamma
gamma = 1.4;

global V1 V2 V3 V4
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

fxS1 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(uL,uR);
fxS2 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(uL,uR).^2 + pavg(rhoL,uL,vL,EL,rhoR,uR,vR,ER);
fxS3 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(uL,uR).*avg(vL,vR);
fxS4 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) (plogmean(rhoL,uL,vL,EL,rhoR,uR,vR,ER)/(gamma-1) + pavg(rhoL,uL,vL,EL,rhoR,uR,vR,ER) + .5*logmean(rhoL,rhoR).*vnormavg(uL,vL,uR,vR)).*avg(uL,uR);

fyS1 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(vL,vR);
fyS2 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(uL,uR).*avg(vL,vR);
fyS3 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) logmean(rhoL,rhoR).*avg(vL,vR).^2 + pavg(rhoL,uL,vL,EL,rhoR,uR,vR,ER);
fyS4 = @(rhoL,uL,vL,EL,rhoR,uR,vR,ER) (plogmean(rhoL,uL,vL,EL,rhoR,uR,vR,ER)/(gamma-1) + pavg(rhoL,uL,vL,EL,rhoR,uR,vR,ER) + .5*logmean(rhoL,rhoR).*vnormavg(uL,vL,uR,vR)).*avg(vL,vR);

% entropy potentials
psix = @(rho,rhou,rhov,E) (gamma-1)*rhou;
psiy = @(rho,rhou,rhov,E) (gamma-1)*rhov;


%% problem params setup

x0 = 0; y0 = 0;

[rhoq uq vq pq] = vortexSolution(xq,yq,0);
rho  = Pq*rhoq;
rhou = Pq*(rhoq.*uq);
rhov = Pq*(rhoq.*vq);
E    = Pq*(pq/(gamma-1) + .5*rhoq.*(uq.^2+vq.^2));

rho = Vq*rho;
rhou = Vq*rhou;
rhov = Vq*rhov;
E = Vq*E;


% xf = reshape(Vfq*x,Nfp,Nfaces*K);
% yf = reshape(Vfq*y,Nfp,Nfaces*K);
% xfnc = reshape(Vfsplit*xf(:,ncfaces),Nfp,num_nonconf_faces*2);
% yfnc = reshape(Vfsplit*yf(:,ncfaces),Nfp,num_nonconf_faces*2);
% eK = find(activeK);
% text(mean(x(:,eK)),mean(y(:,eK)),num2str(eK))
% hold on
% plot([xf xfnc],[yf yfnc],'o')
% % xf(:,ncfaces) = Pqsplit*reshape(xfnc,2*Nfp,num_nonconf_faces);
% % yf(:,ncfaces) = Pqsplit*reshape(yfnc,2*Nfp,num_nonconf_faces);
% % plot(xf,yf,'^')
% return

%%

global wJq
wJq = diag(wq)*(J);

% Runge-Kutta residual storage
res1 = zeros(Nq,K);
res2 = zeros(Nq,K);
res3 = zeros(Nq,K);
res4 = zeros(Nq,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
CNh = max(CN*max(sJ(:)./Jf(:)));
dt = CFL*2/CNh;

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

figure(1)
for i = 1:Nsteps
    for INTRK = 1:5                
        
        % ignore curvilinear for now                                 
        
        % interpolate entropy vars to face points
        q1 = V1(rho,rhou,rhov,E);
        q2 = V2(rho,rhou,rhov,E);
        q3 = V3(rho,rhou,rhov,E);
        q4 = V4(rho,rhou,rhov,E);        
        q1M = VfPq*q1; 
        q2M = VfPq*q2;
        q3M = VfPq*q3;
        q4M = VfPq*q4; 
        
        % evaluate 
        rhoM  = U1(q1M,q2M,q3M,q4M);
        rhouM = U2(q1M,q2M,q3M,q4M);
        rhovM = U3(q1M,q2M,q3M,q4M);
        EM    = U4(q1M,q2M,q3M,q4M);
        
        u = rhou./rho; 
        v = rhov./rho; 
        uM = rhouM./rhoM;
        vM = rhovM./rhoM;
        
        % extra LF flux info        
        [rhs1 rhs2 rhs3 rhs4]  = RHS2Dsimple(rho,u,v,E,rhoM,uM,vM,EM);        
        
        if (INTRK==5)
            rhstest(i) = 0;
            for e = 1:K                
                q1e = J(:,e).*wq.*(q1(:,e));
                q2e = J(:,e).*wq.*(q2(:,e));
                q3e = J(:,e).*wq.*(Vq*q3(:,e));
                q4e = J(:,e).*wq.*(Vq*q4(:,e));
                r1 = rhs1(:,e);
                r2 = rhs2(:,e);
                r3 = rhs3(:,e);
                r4 = rhs4(:,e);
                rhstest(i) = rhstest(i) + (q1e'*r1 + q2e'*r2 + q3e'*r3 + q4e'*r4);                
            end
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
    entropy(i) = sum(sum(wJq.*Sq));
    
    if  (mod(i,5)==0 || i==Nsteps)
        clf
        pp = rho;
        vv = real(Vp*Pq*pp);
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f, N = %d, K1D = %d',dt*i,N,K1D))
        %                 view(3)
        drawnow
    end    
end

[rq2 sq2 wq2] = Cubature2D(Nq+2);
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
err = wJq2.*((rhoq-rhoex).^2 + (rhouq-rhouex).^2 + (rhovq-rhovex).^2 + (Eq-Eex).^2);
L2err = sqrt(sum(err(:)))

dS = abs(entropy-entropy(1));
dS = entropy + max(abs(entropy));
figure(2)
semilogy(dt*(1:Nsteps),dS,'o--');hold on
semilogy(dt*(1:Nsteps),abs(rhstest),'x--');hold on
legend('dS','rhstest')


function [rhs1 rhs2 rhs3 rhs4] = RHS2Dsimple(rhoq,uq,vq,Eq,rhoM,uM,vM,EM)

Globals2D;

global M Vq Pq Lq Lqf Vfqf Vfq VqPq VqLq
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ nrJq nsJq
global mapPq
global fxS1 fyS1 fxS2 fyS2 fxS3 fyS3 fxS4 fyS4
global pfun beta pavg plogmean vnormavg avg
global Drq Dsq VfPq

global DNr DNs

% combine vol/flux vars
rho = [rhoq; rhoM];
u = [uq; uM];
v = [vq; vM];
E = [Eq; EM];
betaq = beta(rhoq,uq,vq,Eq);
betafq = beta(rhoM,uM,vM,EM);
betaN = [betaq;betafq];

rhoP = rhoM(mapPq);
uP = uM(mapPq);
vP = vM(mapPq);
EP = EM(mapPq);

rhoM = reshape(rhoM,Nfp,Nfaces*K);
uM = reshape(uM,Nfp,Nfaces*K);
vM = reshape(vM,Nfp,Nfaces*K);
EM = reshape(EM,Nfp,Nfaces*K);

% Lax-Friedrichs flux
global gamma tau
unorm2 = (uM.^2+vM.^2);
pM = (gamma-1)*(EM - .5*rhoM.*unorm2);
cvel = sqrt(gamma*pM./rhoM);
lam = sqrt(unorm2)+cvel;
LFc = max(lam(mapPq),lam);

dQ1 = rhoP-rhoM;
dQ2 = rhoP.*uP-rhoM.*uM;
dQ3 = rhoP.*vP-rhoM.*vM;
dQ4 = EP-EM; 

% ------------------ 
% do noncon projection stuff here
global V1 V2 V3 V4
for f = 1:Nfaces*K
end
% ------------------ 

Lf1 = tau*reshape(LFc.*dQ1,Nfp*Nfaces,K).*sJ;
Lf2 = tau*reshape(LFc.*dQ2,Nfp*Nfaces,K).*sJ;
Lf3 = tau*reshape(LFc.*dQ3,Nfp*Nfaces,K).*sJ;
Lf4 = tau*reshape(LFc.*dQ4,Nfp*Nfaces,K).*sJ;

% fix this for noncon
fSf1 = nxJ(:).*fxS1(rhoM(:),uM(:),vM(:),EM(:),rhoP(:),uP(:),vP(:),EP(:)) + nyJ(:).*fyS1(rhoM(:),uM(:),vM(:),EM(:),rhoP(:),uP(:),vP(:),EP(:));
fSf2 = nxJ(:).*fxS2(rhoM(:),uM(:),vM(:),EM(:),rhoP(:),uP(:),vP(:),EP(:)) + nyJ(:).*fyS2(rhoM(:),uM(:),vM(:),EM(:),rhoP(:),uP(:),vP(:),EP(:));
fSf3 = nxJ(:).*fxS3(rhoM(:),uM(:),vM(:),EM(:),rhoP(:),uP(:),vP(:),EP(:)) + nyJ(:).*fyS3(rhoM(:),uM(:),vM(:),EM(:),rhoP(:),uP(:),vP(:),EP(:));
fSf4 = nxJ(:).*fxS4(rhoM(:),uM(:),vM(:),EM(:),rhoP(:),uP(:),vP(:),EP(:)) + nyJ(:).*fyS4(rhoM(:),uM(:),vM(:),EM(:),rhoP(:),uP(:),vP(:),EP(:));

fSf1 = reshape(fSf1,Nfp*Nfaces,K) - .25*Lf1;
fSf2 = reshape(fSf2,Nfp*Nfaces,K) - .25*Lf2;
fSf3 = reshape(fSf3,Nfp*Nfaces,K) - .25*Lf3;
fSf4 = reshape(fSf4,Nfp*Nfaces,K) - .25*Lf4;

divF1 = zeros(size(DNr,1),K);
divF2 = zeros(size(DNr,1),K);
divF3 = zeros(size(DNr,1),K);
divF4 = zeros(size(DNr,1),K);
for e = 1:K
    
    [rhox, rhoy] = meshgrid(rho(:,e));
    [ux, uy] = meshgrid(u(:,e));
    [vx, vy] = meshgrid(v(:,e));
    [betax, betay] = meshgrid(betaN(:,e));
    
    % optimized evaluations
    rholog = logmean(rhox,rhoy);
    rhoavg = avg(rhox,rhoy);
    uavg = avg(ux,uy);
    vavg = avg(vx,vy);
    vnavg = 2*(uavg.^2 + vavg.^2) - (avg(ux.^2,uy.^2) + avg(vx.^2,vy.^2));
    pa = rhoavg./(2*avg(betax,betay));
    
    FxS1 = rholog.*uavg;      FyS1 = rholog.*vavg;
    FxS2 = FxS1.*uavg + pa;   FyS2 = FyS1.*uavg;
    FxS3 = FyS2;              FyS3 = FyS1.*vavg + pa;
    f4aux = rholog./(2*(gamma-1)*logmean(betax,betay)) + pa + .5*rholog.*vnavg;
    FxS4 = f4aux.*uavg;
    FyS4 = f4aux.*vavg;
    
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

PN = [VqPq VqLq]; 
rhs1 = 2*PN*divF1 + VqLq*(fSf1);
rhs2 = 2*PN*divF2 + VqLq*(fSf2);
rhs3 = 2*PN*divF3 + VqLq*(fSf3);
rhs4 = 2*PN*divF4 + VqLq*(fSf4);

% collocation assumption
rhs1 = -(rhs1./J);
rhs2 = -(rhs2./J);
rhs3 = -(rhs3./J);
rhs4 = -(rhs4./J);

% rhs1 = -VqPq*(rhs1./J);
% rhs2 = -VqPq*(rhs2./J);
% rhs3 = -VqPq*(rhs3./J);
% rhs4 = -VqPq*(rhs4./J);

global activeK
rhs1(:,find(~activeK)) = 0;
rhs2(:,find(~activeK)) = 0;
rhs3(:,find(~activeK)) = 0;
rhs4(:,find(~activeK)) = 0;


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

% rho = (2 + sin(pi*(x - t)));
% u = ones(size(x));
% v = zeros(size(x));
% p = ones(size(x));

if 0
    % pulse condition
    x0 = 5;
    rho = 2 + (abs(x-x0) < 5);
    u = 0*rho;
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
global qt activeK Vsplit Vfsplit

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
