function acousticCurved

% clear -global *
% clear

Globals2D;

N = 3;
useJprojection = 1;

filename = 'Grid/Other/circA01.neu';
[Nv, VX, VY, K, EToV, BCType] = MeshReaderGambitBC2D(filename);


% This builds the nodal DG stuff
StartUp2D;
BuildBCMaps2D;

nref = 2;
for ref = 1:nref
    Refine2D(ones(size(EToV)));
    StartUp2D;
    BuildBCMaps2D;
end

% Store the domain boundary edges
[k,f] = find(BCType);
curved = sort(unique(k));

% Push all boundary faces to unit cylinder
MakeCylinder2D([k,f], 1, 0, 0);
% PlotMesh2D;return
% Set initial conditions
% First 6 modes of eigenmodes with 6 azimuthal periods
alpha = [2.40482555769577, 5.52007811028631, 8.65372791291101, 11.7915, 14.9309]; % https://math.dartmouth.edu/archive/m23f09/public_html/drum.pdf

% choose radial mode
alpha0 = alpha(2); rad = sqrt(x.^2+y.^2);

p = besselj(0, alpha0*rad);
u = zeros(Np, K); v = zeros(Np, K);

FinalTime = 1;
[u,v,p,time] = acousticsCurved2D(u,v,p,FinalTime);

pex = besselj(0,alpha0*rad)*cos(alpha0*time);
err = p-pex;

% clf;color_line3(x,y,err,err,'.');
% title('Error at time T = 1, N = 4')
cinfo = BuildCurvedOPS2D_opt(3*N,useJprojection);

Mhat = inv(V*V');
L2err2 = 0;
for e = 1:K
    MK = Mhat*J(1,e);
    eK = err(:,e);
    L2err2(e) = eK'*MK*eK;
end
for ee = 1:length(curved)
    MK = cinfo(ee).MM;
    e = curved(ee);
    eK = err(:,e);
    L2err2(e) = eK'*MK*eK;
end
L2err = sqrt(sum(L2err2));
disp(sprintf('L2 error = %e\n',L2err))

% % this estimates computes order of convergence from a vector of err,h
% err = [5.8888e-04, 3.8601e-05, 2.6236e-06];h = [1 .5 .25];
% err = [6.002e-04, 3.89e-05, 2.6838e-06]; % subparam = 2
% fit = [log(h(:)) ones(size(h(:)))]\log(err(:)); order = fit(1)

% keyboard
% L2err1 = [2.429482e-02 4.655834e-03 1.950810e-04 2.093990e-05 8.206020e-07];   % J-dependent mass matrix
% L2err2 = [2.429482e-02 4.655300e-03 1.950115e-04 2.093971e-05 8.203261e-07]; % 1/J-weighted

function [u,v,p,time] = acousticsCurved2D(u, v, p, FinalTime)

% function [u,v,p] = acousticsCurved2D(u, v, p, FinalTime)
% Purpose  : Integrate TM-mode acoustics until FinalTime starting with initial conditions u,v,p

Globals2D;
time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

cinfo = BuildCurvedOPS2D(3*N);

% compute time step size
rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D; dt = min(dtscale)*rmin*2/3

% outer time step loop
tstep = 0;
while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        % compute right hand side of TM-mode acoustics's equations
        [rhsu, rhsv, rhsp] = acousticsCurvedRHS2D(cinfo, u,v,p);
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        p = p+rk4b(INTRK)*resp;
    end;
    
    if 0 && mod(tstep,10)==0
        clf
        color_line3(x,y,p,p,'.');
        view(3)
        axis([-1 1 -1 1 -1.25 1.25])
        drawnow
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
    
end
return

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
        gndotdU = gnx.*gdu + gny.*gdv;
        fluxp = gdp - gndotdU;
        fluxu = (gndotdU - gdp).*gnx;
        fluxv = (gndotdU - gdp).*gny;
        
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
fluxp =  dp - ndotdU;
fluxu =  (ndotdU - dp).*nx;
fluxv =  (ndotdU - dp).*ny;

% local derivatives of fields
[dpdx,dpdy] = Grad2D(p);  [divU] = Div2D(u,v);

% compute right hand sides of the PDE's
rhsp =  -divU + LIFT*(Fscale.*fluxp)/2.0;
rhsu =  -dpdx + LIFT*(Fscale.*fluxu)/2.0;
rhsv =  -dpdy + LIFT*(Fscale.*fluxv)/2.0;

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
curved = sort(unique(k));
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
      
      M = cV'*diag(cW)*cV;
      MinvJ = cV'*diag(cW./cJ)*cV;
      cMM = M*(MinvJ\M);
      
  elseif useJprojection==2 % alias by projecting 1/J onto degree N polynomial first
      
      M = cV'*diag(cW)*cV;
      invcJ = cV*(M\(cV'*(cW./cJ))); % project 1/J onto degree N polynomial, then interp to quadrature           
      MinvJ = cV'*diag(cW.*invcJ)*cV;
      cMM = M*(MinvJ\M);
      
  elseif useJprojection==3 % locally conservative fix
      
      %       disp('using locally conservative version')
      Mex = cV'*diag(cJ.*cW)*cV;
      M = cV'*diag(cW)*cV; % reference M
      
      MinvJ = cV'*diag(cW./cJ)*cV; % quadrature-based 1/J      
%       invcJ = cV*(M\(cV'*(cW./cJ))); % project 1/J onto degree N polynomial, then interp to quadrature
%       MinvJ = cV'*diag(cW.*invcJ)*cV;
      
      
      cMM = M*(MinvJ\M); % non-conservative mass matrix
      eK = ones(Np,1);
      v = (cMM-Mex)*eK; % compute conservation error for each basis fxn
      alpha = 0;
      if abs(v'*eK) > 1e-14
          alpha = 1/(v'*eK);          
      end
      v = sqrt(alpha)*v;      
      cMM = cMM - sign(alpha)*v*v'; % conservative mass matrix      

  else % true mass matrix 
      disp('using true mass')
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