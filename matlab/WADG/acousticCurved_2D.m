function [L2err href] = acousticCurved_2D(Nin,nref,useJprojection)

% clear -global *
% clear

Globals2D;

if nargin==0
    N = 7;
    nref = 2;
    useJprojection = 1;
else
    N = Nin;
end

Nq = 2*N+1;

filename = 'Grid/Other/circA01.neu';
[Nv, VX, VY, K, EToV, BCType] = MeshReaderGambitBC2D(filename);

% K1D = 2;
% [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

% This builds the nodal DG stuff
StartUp2D;

BCType = zeros(size(EToV));
for e = 1:K
    BCType(e,EToE(e,:)==e) = 6;
end

BuildBCMaps2D;
for ref = 1:nref
    Refine2D(ones(size(EToV)));
    StartUp2D;
    BuildBCMaps2D;
end

% Store the domain boundary edges
[k,f] = find(BCType);

% Push all boundary faces to unit cylinder
if (useJprojection >= 0)
    MakeCylinder2D([k,f], 1, 0, 0);    
end
cinfo = BuildCurvedOPS2D_opt(Nq,useJprojection);

% PlotMesh2D

% jesse custom curved driver
global Vq wq Vrq Vsq
global Vfq wfq VfqFace
global rxq sxq ryq syq Jq Jfq nxq nyq sJq
global Pq Pfq


[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

[rq sq wq] = Cubature2D(Nq);
Vq = Vandermonde2D(N,rq,sq)/V;
[Vrq Vsq] = GradVandermonde2D(N,rq,sq);
Vrq = Vrq/V; Vsq = Vsq/V;

% interp nodal to quadrature
[rq1D, wq1D] = JacobiGQ(0, 0, Nq);
rfq = [rq1D, -rq1D, -ones(size(rq1D))];
sfq = [-ones(size(rq1D)), rq1D, -rq1D];
rfq = rfq(:); sfq = sfq(:);
wfq = repmat(wq1D,Nfaces,1);
Vfq = Vandermonde2D(N,rfq,sfq)/V;
[Vrfq Vsfq] = GradVandermonde2D(N,rfq,sfq);
Vrfq = Vrfq/V;
Vsfq = Vsfq/V;

V1D = Vandermonde1D(N,JacobiGL(0,0,N));
VfqFace = Vandermonde1D(N,rq1D)/V1D;
VfqFaceperm = Vandermonde1D(N,rq1D(end:-1:1))/V1D;
%VfqFace = blkdiag(VfqFace,VfqFace,VfqFaceperm); % repeat for 3 faces
VfqFace = blkdiag(VfqFace,VfqFace,VfqFace); % repeat for 3 faces


% project down
Pq = (V*V')*(Vq'*diag(wq));
Pfq = (V*V')*(Vfq'*diag(wfq));

Nc = length(wq);
rxq = zeros(Nc,K); sxq = zeros(Nc,K);
ryq = zeros(Nc,K); syq = zeros(Nc,K);
Jq = zeros(Nc,K);
Nfc = length(cinfo(1).gnx(:));
nxq = zeros(Nfc,K); nyq = zeros(Nfc,K);
sJq = zeros(Nfc,K);
Jfq = zeros(Nfc,K);
for e = 1:K
    [rxk,sxk,ryk,syk,Jk] = GeometricFactors2D(x(:,e),y(:,e),Vrq,Vsq);
    rxq(:,e) = rxk;
    sxq(:,e) = sxk;
    ryq(:,e) = ryk;
    syq(:,e) = syk;
    Jq(:,e) = Jk;
    
    nxq(:,e) = cinfo(e).gnx(:);
    nyq(:,e) = cinfo(e).gny(:);
    sJq(:,e) = cinfo(e).gsJ(:);
    [~,~,~,~,JfqK] = GeometricFactors2D(x(:,e),y(:,e),Vrfq,Vsfq);
Jfq(:,e) = JfqK;
end


global tau;
tau = 1;


%%
if 0 && nargin==0 
    
    e = zeros(3*Np*K,1);
    A = zeros(3*Np*K);
    for i = 1:3*Np*K
        e(i) = 1;
        ids = 1:Np*K;
        p = reshape(e(ids),Np,K);
        u = reshape(e(ids + Np*K),Np,K);
        v = reshape(e(ids + 2*Np*K),Np,K);
        if useJprojection==1
            [rhsu, rhsv, rhsp] = acousticsRHS2D_JC( u,v,p);
        elseif useJprojection==0
            [rhsu, rhsv, rhsp] = acousticsCurvedRHS2D(cinfo, u,v,p);
        else
            [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p); 
        end
        
        A(:,i) = [rhsp(:);rhsu(:);rhsv(:)];
        e(i) = 0;
        if mod(i,100)==0
            disp(sprintf('on column %d out of %d\n',i,3*Np*K))
        end
    end
    
%     [rq sq wq] = Cubature2D(Nq);
%     invMhat = V*V';
%     for e = 1:K                 
%         invMK = invMhat*Vq'*diag(wq./Jq(:,e))*Vq*invMhat;
%         Mblk{e} = inv(invMK); %kron(eye(3),MK);
%     end
%     Mblk = blkdiag(Mblk{:});
%     Mblk = blkdiag(Mblk,Mblk,Mblk); % for 3 fields
%     MA = Mblk*A;
    
    lam = eig(A);
    hold on;
    if useJprojection==1
        plot(lam,'o')
    else
        plot(lam,'*')
    end
    title(sprintf('Largest real part = %e',max(real(lam))))
    return

end
%%

% Set initial conditions
% First 6 modes of eigenmodes with 6 azimuthal periods
alpha = [2.40482555769577, 5.52007811028631, 8.65372791291101, 11.7915, 14.9309]; % https://math.dartmouth.edu/archive/m23f09/public_html/drum.pdf

% choose radial mode
alpha0 = alpha(2); rad = sqrt(x.^2+y.^2);

p = besselj(0, alpha0*rad);

d = 25; x0 = -1/3; y0 = 1/3;
p = exp(-d*sqrt((x-x0).^2 + (y-y0).^2)) + exp(-10*sqrt((x+x0).^2 + (y+y0).^2));

u = zeros(Np, K); v = zeros(Np, K);

FinalTime = 2;

% setup
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);
rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D; 

dt = .75*min(dtscale)*rmin;

% outer time step loop
time = 0; tstep = 0;
while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        % compute right hand side of TM-mode acoustics's equations
%         
        if (useJprojection==0)
            [rhsu, rhsv, rhsp] = acousticsCurvedRHS2D(cinfo, u,v,p);            
        elseif (useJprojection==-1)
            [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p);        
        else
            [rhsu, rhsv, rhsp] = acousticsRHS2D_JC(u,v,p);        
        end
%         
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        p = p+rk4b(INTRK)*resp;
    end;
    
    if 1 && mod(tstep,10)==0
        clf
        pp = Vp*p;
        color_line3(xp,yp,pp,pp,'.');
        view(2)
%         axis([-1 1 -1 1 -.5 .5])
        title(sprintf('time = %f',time));
        drawnow
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
    if mod(tstep,10) ==0
        disp(sprintf('on timestep %d out of %d\n',tstep,round(FinalTime/dt)))
    end
end

xq = Vq*x; 
yq = Vq*y;
radq = sqrt(xq.^2 + yq.^2);
pex = besselj(0,alpha0*radq)*cos(alpha0*time);

L2err = 0.0;
for e = 1:K
    errK = sum(wq.*Jq(:,e).*(Vq*p(:,e)-pex(:,e)).^2);
    L2err = L2err + errK;
end
L2err = sqrt(L2err);
href = max(1./Fscale(:));
disp(sprintf('L2 error = %e\n',L2err))


function [rhsu, rhsv, rhsp] = acousticsCurvedRHS2D(cinfo, u,v,p)

% function [rhsu, rhsv, rhsp] = acousticsCurvedRHS2D(cinfo, u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

%SubParametric for all non-curved
[rhsu,rhsv,rhsp] = acousticsRHS2D(u, v, p);

Ncinfo = length(cinfo);

% correct residuals at each curved element
for n=1:Ncinfo
    
    % for each curved element computed L2 derivatives via cubature
    cur = cinfo(n); k1 = cur.elmt; cDx = cur.Dx; cDy = cur.Dy;
    
    rhsp(:,k1) =  -(cDx*u(:,k1) + cDy*v(:,k1)); % -div U
    rhsu(:,k1) =  -cDx*p(:,k1);  % -grad p
    rhsv(:,k1) =  -cDy*p(:,k1);
    
    % for each face of each curved element use Gauss quadrature based lifts
    for f1=1:Nfaces
        k2 = EToE(k1,f1);
        gnx = cur.gnx(:,f1); gny = cur.gny(:,f1);
        gVM = cur.gVM(:,:,f1); gVP = cur.gVP(:,:,f1);
        glift = cur.glift(:,:,f1);
        
        gdp = gVP*p(:,k2) - gVM*p(:,k1); % [p] = p+ - p-
        gdu = gVP*u(:,k2) - gVM*u(:,k1);
        gdv = gVP*v(:,k2) - gVM*v(:,k1);
        
        % compute difference of solution traces at Gauss nodes
        % correct jump at Gauss nodes on domain boundary faces
        if(k1==k2)
            gdu = 0*gdu; gdv = 0*gdv;
            gdp = -2.0*gVM*p(:,k1); %p+ = -p-
        end
        
        % perform upwinding
        global tau
        gndotdU = gnx.*gdu + gny.*gdv;
        fluxp = tau*gdp - gndotdU;
        fluxu = (tau*gndotdU - gdp).*gnx;
        fluxv = (tau*gndotdU - gdp).*gny;
        
        % lift flux terms using Gauss based lift operator
        rhsp(:,k1) = rhsp(:,k1) + glift*fluxp/2;
        rhsu(:,k1) = rhsu(:,k1) + glift*fluxu/2;
        rhsv(:,k1) = rhsv(:,k1) + glift*fluxv/2;
    end
end


return;

function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

% Impose reflective boundary conditions (p+ = -p-)
du(mapB) = 0; dv(mapB) = 0; dp(mapB) = -2*p(vmapB);

% evaluate upwind fluxes
ndotdU = nx.*du + ny.*dv;
global tau
fluxp =  tau*dp - ndotdU;
fluxu =  (tau*ndotdU - dp).*nx;
fluxv =  (tau*ndotdU - dp).*ny;

% local derivatives of fields
pr = Dr*p; ps = Ds*p;
ur = Dr*u; us = Ds*u;
vr = Dr*v; vs = Ds*v;

dpdx = pr.*rx + ps.*sx;
dpdy = pr.*ry + ps.*sy;
dudx = ur.*rx + us.*sx;
dvdy = vr.*ry + vs.*sy;
divU = dudx + dvdy; 

% compute right hand sides of the PDE's
rhsp =  -divU + LIFT*(Fscale.*fluxp)/2.0;
rhsu =  -dpdx + LIFT*(Fscale.*fluxu)/2.0;
rhsv =  -dpdy + LIFT*(Fscale.*fluxv)/2.0;

return;


function [rhsu, rhsv, rhsp] = acousticsRHS2D_JC(u,v,p)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

global Vq Vrq Vsq VfqFace
global rxq sxq ryq syq Jq Jfq nxq nyq sJq
global Pq Pfq

Globals2D;

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

% Impose reflective boundary conditions (p+ = -p-)
du(mapB) = 0; dv(mapB) = 0; dp(mapB) = -2*p(vmapB);

% can interp numerical fluxes *after* since mult by sJ is higher order
dp = VfqFace*dp;
du = VfqFace*du;
dv = VfqFace*dv;

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

% compute right hand sides of the PDE's
rhsp =  Pq*(-divU.*Jq) + Pfq*(fluxp.*sJq/2.0);
rhsu =  Pq*(-dpdx.*Jq) + Pfq*(fluxu.*sJq/2.0);
rhsv =  Pq*(-dpdy.*Jq) + Pfq*(fluxv.*sJq/2.0);

% apply inverse mass matrix
rhsp = Pq*((Vq*rhsp)./Jq);
rhsu = Pq*((Vq*rhsu)./Jq);
rhsv = Pq*((Vq*rhsv)./Jq);

return;


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
        
    elseif useJprojection==2 % alias by projecting 1/J onto degree N polynomial first
        
        M = cV'*diag(cW)*cV;
        invcJ = cV*(M\(cV'*(cW./cJ))); % project 1/J onto degree N polynomial, then interp to quadrature
        MinvJ = cV'*diag(cW.*invcJ)*cV;
        cMM = M*(MinvJ\M);
    
    elseif useJprojection==3 % just use reference matrix - should suck.
        
        cMM = cV'*diag(cW)*cV;
        
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