function [cinfo] = BuildCurvedOPS2D_opt(intN,intNf)

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
[gz, gw] = JacobiGQ(0, 0, intNf);

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
    
    % build individual lift matrices at each face
    for f1=1:Nfaces
        
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
        

    end
end

return