clear
hybridgGlobals3D
hybridgGlobalFlags
useSEM = 1; useLSC = 1; useNodalTets = 1; useSkewRHS = 1;

% hybrid_mesh
% prism2; EToV = EToV'; EToE = EToE'; EToF = EToF'; 
% prism_Twist2;
% pyr_mesh6; EToV = EToV'; EToE = EToE'; EToF = EToF'; 
% tet_mesh; 
cube_mesh; 

% VX = [    -1     1    1     -1    -1     1    1     -1]'; 
% VY = [    -1    -1     1     1    -1    -1     1     1]';  
% VZ = [    -1    -1    -1    -1     1     1     1     1]';
% EToV = 1:length(VX); EToE = 0*EToV; EToF = 0*EToV;

% hex_mesh
K = size(EToV,1);
N = 3;

hybridgStartUp

%% set initial conditions

% exact sol
uexf = @(x,y,z,time) cos(pi/2*x).*cos(pi/2.*y).*cos(pi/2.*z)*cos(sqrt(3)*pi/2*time);
f = @(x,y,z) uexf(x,y,z,0);

b = zeros(NpMax,K);
b(1:NpH,hexK ) = VH'*(wJ(1:NcH,hexK ).*f(xqH,yqH,zqH));
b(1:NpW,wedgK) = VW'*(sqJW.*wJ(1:NcW,wedgK).*f(xqW,yqW,zqW));
b(1:NpP,pyrK ) = VP'*(wJ(1:NcP,pyrK ).*f(xqP,yqP,zqP));

p = invM.*b; % LSC inverse built into invM - extra storage but OK

p(1:NpT,tetK ) = f(xT,yT,zT); % interpolate for tets

%% fine spaced plotting - things to load

Nplot = 2;
rp1D = linspace(-1,1,Nplot+1);

% hex
%[rpH spH tpH] = meshgrid(rp1D); 
%rpH = rpH(:); spH = spH(:); tpH = tpH(:);
sk=1;
for i = 1:Nplot+1
    for j = 1:Nplot+1
        for k = 1:Nplot+1
            rpH(sk) = rp1D(k);
            spH(sk) = rp1D(j);
            tpH(sk) = rp1D(i);
            sk = sk + 1;
        end
    end
end

% wedge 
[rpTri spTri] = EquiNodes2D(Nplot); [rpTri spTri] = xytors(rpTri,spTri);
rpW = []; spW = []; tpW = [];
for i = 1:Nplot+1
    rpW = [rpW(:); rpTri(:)];    tpW = [tpW(:); spTri(:)];
    spW = [spW(:); ones(size(rpTri))*rp1D(i)];
end

% pyr
apP = []; bpP = []; cpP = [];
c = linspace(-1,1,Nplot+1);
for level = 0:Nplot
    if level < Nplot        
        a1D = linspace(-1,1,Nplot+1-level);
    else
        a1D = 0;
    end    
    [aa bb] = meshgrid(a1D);
    apP = [apP; aa(:)];  bpP = [bpP; bb(:)];    
    cpP = [cpP; c(level+1)*ones(size(aa(:)))];
end
[rpP spP tpP] = pyr_abctorst(apP,bpP,cpP);

% tet
[rpT spT tpT] = EquiNodes3D(Nplot);

% reference bases
VDMH = hex_basis(N,rpH,spH,tpH);
VDMW = wedge_basis(N,rpW,spW,tpW);
VDMP = pyr_basis(N,rpP,spP,tpP);
VDMT = tet_basis(N,rpT,spT,tpT)/VTnodal;

% reference triangulations index from 0
TH = delaunayFixVol(rpH,spH,tpH)-1;
TW = delaunayFixVol(rpW,spW,tpW)-1;
TP = delaunayFixVol(rpP,spP,tpP)-1;
TT = delaunayFixVol(rpT,spT,tpT)-1;

% hold on
% for t = 1:size(TP,1)
%     tet = TP(t,:)+1; edges = tet([1 2 2 4 4 1 1 3 3 4 4 2 2 3]);
%     plot3(rpP(edges),spP(edges),tpP(edges),'.-')
% end
% return
% 
% hold on
% for t = 1:size(TT,1)
%     tet = TT(t,:)+1; edges = tet([1 2 2 4 4 1 1 3 3 4 4 2 2 3]);
%     plot3(rpT(edges),spT(edges),tpT(edges),'.-')
% end
% return
%% physical points

xp = zeros(NpMax,K); yp = zeros(NpMax,K); zp = zeros(NpMax,K); pplot = zeros(NpMax,K);

for e = 1:K
    v = EToV(e,:); v = v(v > 0); NvK = nnz(v);
    
    switch NvK
        case 4 % tet
            [xpK,ypK,zpK] = tet_geom_factors(VX(v),VY(v),VZ(v),rpT,spT,tpT);
            pplotK = VDMT * p(1:NpT,e);
            pplot(1:length(rpT),e) = pplotK;
            NpK = length(rpT);
            
        case 5 % pyr
            
            [xpK,ypK,zpK] = pyr_geom_factors(VX(v),VY(v),VZ(v),rpP,spP,tpP);
            pplot(1:length(rpP),e) = VDMP * p(1:NpP,e);
            NpK = length(rpP);
            
        case 6 % wedge
            
            [xpK,ypK,zpK, ~,~,~,~,~,~,~,~,~,JpW] = ...
                wedge_geom_factors(VX(v),VY(v),VZ(v),rpW,spW,tpW);
            pplotK = (VDMW * p(1:NpW,e))./sqrt(JpW);
            pplot(1:length(rpW),e) = pplotK;
            NpK = length(rpW);
            
        case 8 % hex
            [xpK,ypK,zpK] = hex_geom_factors(VX(v),VY(v),VZ(v),rpH,spH,tpH);
            pplotK = VDMH * p(1:NpH,e);
            pplot(1:length(rpH),e) = pplotK;
            NpK = length(rpH);
            
    end
    xp(1:NpK,e) = xpK;    yp(1:NpK,e) = ypK;    zp(1:NpK,e) = zpK;        
end
color_line3(xp,yp,zp,pplot,'.')

%% write to paraview
fID = fopen(sprintf('p_output.vtk'),'w');

fprintf(fID,'# vtk DataFile Version 3.0\n');
fprintf(fID,'hybridg paraview-compatible vtk output\n');
fprintf(fID,'ASCII\n');
fprintf(fID,'DATASET UNSTRUCTURED_GRID\n');

% xyz points total
Npts = length(hexK)*length(rpH) + length(wedgK)*length(rpW) + ...
    length(pyrK)*length(rpP) + length(tetK)*length(rpT);

fprintf(fID,'POINTS %d float\n',Npts);
NplotH = length(rpH); 
for ee = 1:length(hexK)
    e = hexK(ee);
    for i = 1:NplotH
        fprintf(fID,'%f %f %f\n',xp(i,e), yp(i,e), zp(i,e));
    end
end

NplotW = length(rpW);
for ee = 1:length(wedgK)
    e = wedgK(ee);
    for i = 1:NplotW
        fprintf(fID,'%f %f %f\n',xp(i,e), yp(i,e), zp(i,e));
    end
end

NplotP = length(rpP);
for ee = 1:length(pyrK)
    e = pyrK(ee);
    for i = 1:NplotP
        fprintf(fID,'%f %f %f\n',xp(i,e), yp(i,e), zp(i,e));
    end
end

NplotT = length(rpT);
for ee = 1:length(tetK)
    e = tetK(ee);
    for i = 1:NplotT
        fprintf(fID,'%f %f %f\n',xp(i,e), yp(i,e), zp(i,e));
    end
end
fprintf(fID,'\n');


% triangulation of each element 
NcellsH = size(TH,1) * length(hexK);
NcellsW = size(TW,1) * length(wedgK);
NcellsP = size(TP,1) * length(pyrK);
NcellsT = size(TT,1) * length(tetK);
Ncells = NcellsH + NcellsW + NcellsP + NcellsT;

fprintf(fID,'CELLS %d %d\n',Ncells, Ncells*5); % 5 pieces of data for a tet
for ee = 1:length(hexK)
    TK = TH + (ee-1)*NplotH; % index by # cells
    for i = 1:size(TH,1)
        fprintf(fID,'4 %d %d %d %d\n', TK(i,:)); % index cells
    end
end

nptsH = length(hexK)*NplotH;
nptsW = length(wedgK)*NplotW;
nptsP = length(pyrK)*NplotP;

offset = nptsH;
for ee = 1:length(wedgK)
    TK = TW + (ee-1)*NplotW; % index by # cells
    for i = 1:size(TW,1)
        fprintf(fID,'4 %d %d %d %d\n', TK(i,:) + offset); % index cells
    end
end

offset = nptsH + nptsW;
for ee = 1:length(pyrK)
    TK = TP + (ee-1)*NplotP; % index by # cells
    for i = 1:size(TP,1)
        fprintf(fID,'4 %d %d %d %d\n', TK(i,:) + offset); % index cells
    end
end

offset = nptsH + nptsW + nptsP;
for ee = 1:length(tetK)
    TK = TT + (ee-1)*NplotT; % index by # cells
    for i = 1:size(TT,1)
        fprintf(fID,'4 %d %d %d %d\n', TK(i,:) + offset); % index cells
    end
end
fprintf(fID,'\n');

fprintf(fID,'CELL_TYPES %d\n',Ncells); 
fprintf(fID,'%d\n', 10*ones(Ncells,1)); % all cells = tets
fprintf(fID,'\n');

% 10 = linear tet: we just tetrahedralize everything.
fprintf(fID,'POINT_DATA %d\n',Npts); 
fprintf(fID,'SCALARS var_name float\n',Npts); 
fprintf(fID,'LOOKUP_TABLE default\n');

for ee = 1:length(hexK)
    e = hexK(ee);
    fprintf(fID,'%f\n',pplot(1:NplotH,e));    
end
for ee = 1:length(wedgK)
    e = wedgK(ee);
    fprintf(fID,'%f\n',pplot(1:NplotW,e));    
end
for ee = 1:length(pyrK)
    e = pyrK(ee);
    fprintf(fID,'%f\n',pplot(1:NplotP,e));    
end
for ee = 1:length(tetK)
    e = tetK(ee);
    fprintf(fID,'%f\n',pplot(1:NplotT,e));    
end
fprintf(fID,'\n');
fclose(fID);

