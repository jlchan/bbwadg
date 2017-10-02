%clear
function [VX VY VZ EToV EToE EToF K num_triangles] = smoothed_wedge_layers(filename,surface_mesh,usurf,usurf2,nlevel)

if nargin==0
    filename = 'cube';
end
if nargin<2
    surface_mesh = '../Grid/Maxwell2D/Maxwell00625.neu';
end
if nargin<3
    %     usurf = @(x,y) .25*sin(2*x+2.5).*sin(3*y+3.5);
    usurf2 = @(x,y) ones(size(x)); % assumes wavy layer is centered around z = 0.
    %     usurf = @(x,y) .5*(sin(2*x+2.5).*sin(3*y+3.5) + .25*sin(4*x+1.5).*sin(5*y+5.5) + .125*sin(8*x).*sin(8*y));
    %usurf = @(x,y) .25*(sin(2*x+2.5).*sin(3*y+3.5) + .25*sin(4*x+1.5).*sin(5*y+5.5) + .125*sin(8*x).*sin(8*y));
    %usurf2 = @(x,y) 1+.1*sin(x+2.5).*sin(3*y+3.5) + .05*cos(10*x).*sin(9*y);
    usurf = @(x,y) -ones(size(x)); % assumes wavy layer is centered around z = 0.
    
end

if nargin<4
    %     nlevel = 1;
    nlevel = 1;
end

% smoothing of triangular surface
Globals2D

% Polynomial order used for approximation
N = 1;

% Read in Mesh
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D(surface_mesh);

num_triangles = K;

% % plotting mesh`
% xe = repmat(x,2,1); ye = repmat(y,2,1);
% ids = [1 2 3 1 4 5 6 4 5 2 3 6];

% Initialize solver and construct grid and metric
StartUp2D;

[R vmapBT] = getCGRestriction();

u = usurf(x,y);
tri = delaunay(x,y);

% init 3D mesh objects
VX3D = VX; VY3D = VY;
VZ3D = diag(1./sum(R,2))*R*u(:); 
VZ3D = VZ3D(:)';
% VX3D = []; VY3D = []; VZ3D = [];

EToV3D = [];
EToF3D = zeros(K*nlevel,5);
EToE3D = zeros(K*nlevel,5);

% points to evaluate at
% yy = linspace(0,1,nlevel+1);
%  yy = linspace(-1,1,nlevel+1); % for wedges only

nlevel
[r1 s1 t1] = wedge_nodes(1);
t1D = linspace(-1,1,nlevel+1);

Koff = 0; Nvoff = Nv;
for level = 1:nlevel
    
    if level==1
        VX2D = VX(:)';
        VY2D = VY(:)';
        uVZ = diag(1./sum(R,2))*R*u(:); % global vertex values
        VZ2D = uVZ(:)';      
    end
    
    
    % plotting
    if nargin==0
        hold on; trisurf(tri,x,y,u)
    end
    
    if 0
        % simple decay map
        f = @(y) (1-y); % seems to work well.  can change exponent
        u0 = u(:) - mean(u(:));
        u(:) = yy(level+1) + u0(:)*f(yy(level+1));
    else % interpolation
        [r s] = Nodes2D(1); [r s] = xytors(r,s);
        [r t] = meshgrid(r,t1D(level+1)); [s, ~] = meshgrid(s,t1D(level+1));
        r = r(:); s = s(:); t = t(:);
        
        %plot3(r,s,t,'o');hold on;
        VW = wedge_basis(1,r,t,s)/wedge_basis(1,r1,s1,t1);
        u = VW*[usurf(x,y); usurf2(x,y)];
        
    end    
    
    % boilerplate stuff
    VX3D = [VX3D VX];
    VY3D = [VY3D VY];
    uVZ = diag(1./sum(R,2))*R*u(:); %u(:);
    VZ3D = [VZ3D uVZ(:)'];
    
    % layer of elements above = Nv verts up
    newEToV = [EToV+Nvoff-Nv (EToV+Nvoff)];
    EToV3D = [EToV3D; newEToV];
    
    % build connectivity based on 2D connectivity
    eids = Koff + (1:K);
    EToE3D(eids,2:4) = EToE + Koff;
    EToE3D(eids,1) = eids - K; % lower layer
    EToE3D(eids,5) = eids + K; % upper layer
    if level==1
        EToE3D(eids,1) = eids; % elem is attached to itself
    end
    if level == nlevel
        EToE3D(eids,5) = eids; % elem is attached to itself
    end
    
    EToF3D(eids,2:4) = EToF + 1; % offset since 1st face is now bottom
    EToF3D(eids,[1 5]) = repmat([5 1],K,1);
    
    % move to next layer of elements
    Koff = Koff + K; % keeps track of how many elements added
    Nvoff = Nvoff + Nv; % number of vertices
    
end

% plotting
if nargin==0
    hold on; trisurf(tri,x,y,u)
end
% ue = [u;uold]; plot3(xe(ids,:),ye(ids,:),ue(ids,:),'ko-','linewidth',2);hold on
% trisurf(tri,x,y,u)
% axis equal

if nargin==0
    return
end

num_surface_verts = length(VX);
EToV2D = EToV;

% convert to 3D mesh quantities
VX = VX3D; VY = VY3D; VZ = VZ3D;
EToV = EToV3D;
EToE = EToE3D;
EToF = EToF3D;
K = size(EToV,1);
save([filename '.mat'],'VX','VY','VZ','EToV','EToE','EToF','K')

%% write out surface triangulation to tetgen file

% add bottom square vertices
VXtg = [-1 1 1 -1 VX2D];
VYtg = [-1 -1 1 1 VY2D];
VZtg = [-1 -1 -1 -1 VZ2D];
EToV2D = EToV2D + 4; % increment all counts by 4 (extra nodes added)

tetgenfile = ['Grid/' filename '.poly'];
fid = fopen(sprintf(tetgenfile), 'w');
fprintf(fid, '# Part 1: node list\n');
fprintf(fid, sprintf('%d 3 0 0\n', 4 + num_surface_verts)); % number of nodes, # dimensions, #attributes, boundary marker 0 or 1
fprintf(fid,'# Node index, node coordinates\n');
VXYZ = [VXtg(:) VYtg(:) VZtg(:)];
% node_ids = 1:size(VXYZ,1);
% fprintf(fid, '%d %14.25f %14.25f %14.25f\n',node_ids(:),VXYZ(:,1:3))
for i = 1:size(VXYZ,1)
    fprintf(fid, '%d %f %f %f\n',i,VXYZ(i,1:3));
end

fprintf(fid, sprintf('\n# Part 2: facet list\n'));
% number of facets, # of holes in facet, boundary marker

fprintf(fid,'# facet count, no boundary marker\n');
fprintf(fid,'%d 0\n',5+size(EToV2D,1));
fprintf(fid,'# facets\n');

% face 1: t = -1
fprintf(fid, '1\n');  % bottom square
fprintf(fid, '4 1 2 3 4\n');  % bottom square
face = [1 2 3 4];
% plot3(VXtg(face),VYtg(face),VZtg(face),'o');hold on
% text(VXtg(face),VYtg(face),VZtg(face),num2str((1:length(face))'))

% face 2: s = -1
fids = find(abs(VYtg+1)<1e-8 & VZtg > -1);
[~,p] = sort(VXtg(fids),'descend'); fids = fids(p);
face = [1 2 fids];
fprintf(fid, '1\n');
fprintf(fid, '%d ',length(face(:)));
for i = 1:length(face(:))
    fprintf(fid, '%d ',face(i));
end
fprintf(fid,'\n');
% plot3(VXtg(face),VYtg(face),VZtg(face),'o');hold on
% text(VXtg(face),VYtg(face),VZtg(face),num2str((1:length(face))'))

% face 3: r = 1
fids = find(abs(VXtg-1)<1e-8 & VZtg > -1);
[~,p] = sort(VYtg(fids),'descend'); fids = fids(p);
face = [2 3 fids];
fprintf(fid, '1\n');
fprintf(fid, '%d ',length(face(:)));
for i = 1:length(face(:))
    fprintf(fid, '%d ',face(i));
end
fprintf(fid,'\n');
% plot3(VXtg(face),VYtg(face),VZtg(face),'ro');hold on
% text(VXtg(face),VYtg(face),VZtg(face),num2str((1:length(face))'))

% face 4: s = 1
fids = find(abs(VYtg-1)<1e-8 & VZtg > -1);
[~,p] = sort(VXtg(fids),'ascend'); fids = fids(p);
face = [3 4 fids];
fprintf(fid, '1\n');
fprintf(fid, '%d ',length(face(:)));
for i = 1:length(face(:))
    fprintf(fid, '%d ',face(i));
end
fprintf(fid,'\n');
% plot3(VXtg(face),VYtg(face),VZtg(face),'mo');hold on
% text(VXtg(face),VYtg(face),VZtg(face),num2str((1:length(face))'))

% face 5: r = -1
fids = find(abs(VXtg+1)<1e-8 & VZtg > -1);
[~,p] = sort(VYtg(fids),'ascend'); fids = fids(p);
face = [4 1 fids];
fprintf(fid, '1\n');
fprintf(fid, '%d ',length(face(:)));
for i = 1:length(face(:))
    fprintf(fid, '%d ',face(i));
end
fprintf(fid,'\n');
% plot3(VXtg(face),VYtg(face),VZtg(face),'mo');hold on
% text(VXtg(face),VYtg(face),VZtg(face),num2str((1:length(face))'))

% treat each surface triangle as a face
for e = 1:size(EToV2D,1)
    
    face = EToV2D(e,:);
    fprintf(fid, '1\n');
    fprintf(fid, '%d ',length(face(:)));
    for i = 1:length(face(:))
        fprintf(fid, '%d ',face(i));
    end
    fprintf(fid,'\n');
    
    %     plot3(VXtg(face),VYtg(face),VZtg(face),'mo');hold on
    %     text(VXtg(face),VYtg(face),VZtg(face),num2str(face(:)))
    %     hold on;
    %     title(sprintf('Face %d',e))
    %     pause
    %     clf
    
end

fprintf(fid, sprintf('\n# Part 3: hole list\n'));
fprintf(fid, sprintf('0\n'));
fprintf(fid, sprintf('\n# Part 4: region list\n'));
fprintf(fid, sprintf('0\n'));

% fprintf(fid, sprintf('\n#endif\n'));
fclose(fid);

% plot3(VXtg,VYtg,VZtg,'o')
% return
%
% drawWedgeMesh(VX,VY,VZ,EToV)


function [Nv, VX, VY, K, EToV] = MeshReaderGmsh2D(FileName)

% function [Nv, VX, VY, K, EToV] = MeshReaderGmsh2D(FileName)
% Purpose  : Read in basic grid information to build grid

file = fopen(FileName, 'rt');
line = fgetl(file);

while(line ~= -1)
    
    if strfind(line,'$Nodes')
        
        line = fgetl(file);
        Nv = str2num(line);
        VX = zeros(1,Nv);
        VY = zeros(1,Nv);
        
        for i=1:Nv
            line = fgetl(file);
            out = sscanf(line,'%i%f%f%f');
            VX(i) = out(2);
            VY(i) = out(3);
        end
    end
    
    if strfind(line,'$Elements')
        
        line = fgetl(file);
        K = str2num(line);
        EToV = zeros(K,3);
        
        kTri = 0;
        for i=1:K
            line = fgetl(file);
            out = sscanf(line,'%i%i');
            if out(2) == 2
                kTri = kTri+1;
                out = sscanf(line,'%i%i%i%i%i%i%i%i');
                EToV(kTri,1) = out(6);
                EToV(kTri,2) = out(7);
                EToV(kTri,3) = out(8);
            end
        end
        
        K = kTri;
        EToV = EToV(1:K,1:3);
    end
    
    line = fgetl(file);
end

