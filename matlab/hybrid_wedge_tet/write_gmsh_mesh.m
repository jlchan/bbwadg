% filename='gmsh_simple';
filename='wavy';
% filename = 'wavy_layer_ref';
ref = 6;
surface_mesh = sprintf('Grid/tri%d.neu',ref);
usurf = @(x,y) .25*(sin(2*x+2.5).*sin(3*y+3.5) + .25*sin(4*x+1.5).*sin(5*y+5.5) + .125*sin(8*x).*sin(8*y));
usurf2 = @(x,y) 1+.1*sin(x+2.5).*sin(3*y+3.5) + .05*cos(10*x).*sin(9*y);
level = [1 1 2.^(1:5)];
nlevel = 1; level(ref);
[VXW VYW VZW EToVW EToEW EToFW KW num_interface_elems] = smoothed_wedge_layers(filename,surface_mesh,usurf,usurf2,nlevel);

% run tetgen to get mesh
file = ['Grid/' filename];
h = [100 1 .01 .005 .001 .0001]; % heuristic! 
maxvol = h(ref);
system(['"./iso2mesh/bin/tetgen.mexmaci64" -pq1.414Va' num2str(maxvol) ' "' file '"'])
% system(['"/Users/jlchan/Desktop/total_hybrid/matlab/TOTAL_wedge/iso2mesh/bin/tetgen.mexmaci64" ' num2str(maxvol) ' "' file '"'])

% read in the generated mesh
[node,elem,face] = readtetgen([file '.1']);

VXT = node(:,1); VYT = node(:,2); VZT = node(:,3);
EToVT = elem;

% plotmesh(node(:,5:3),elem);
% % drawTetMesh(VXT,VYT,VZT,EToVT);hold on
% drawTetWedgeMesh(VXW,VYW,VZW,EToVW);axis equal
% return

%% 
EToV = [EToVT];

% % combine nodes, trim one copy of interface nodes
VX = [VXT(:)'];
VY = [VYT(:)'];
VZ = [VZT(:)'];

write_gmsh_file(filename,VX,VY,VZ,EToV)

return

%% combine list of vertices

num_interface_nodes = length(VXW)/(nlevel+1);
VXI = VXW(1:num_interface_nodes); % first nodes = interface nodes
VYI = VYW(1:num_interface_nodes);
VZI = VZW(1:num_interface_nodes); 

% [~,p] = sort(VZT,'descend');
% interface_ids = p(1:num_interface_nodes);
interface_ids = (1:num_interface_nodes)+4; % from adding 4 extra box nodes to the tet mesher

% plot3(VXT(interface_ids),VYT(interface_ids),VZT(interface_ids),'bo');hold on
% plot3(VXW(1:num_interface_nodes),VYW(1:num_interface_nodes),VZW(1:num_interface_nodes),'r*')
% text(VXW(1:num_interface_nodes)+.1,VYW(1:num_interface_nodes),VZW(1:num_interface_nodes),num2str((1:length(VXW(1:num_interface_nodes)))'))
% return

% find matches in interface nodes
tmp = ones(1,num_interface_nodes);
xM = VXT(interface_ids)*tmp; yM = VYT(interface_ids)*tmp; zM = VZT(interface_ids)*tmp;
xP = VXI(:)*tmp; yP = VYI(:)*tmp; zP = VZI(:)*tmp;
D = (xM -xP').^2 + (yM-yP').^2 + (zM-zP').^2; 
[i, j] = find(abs(D)<1e-8); 
% (VXT(interface_ids(i)) - VXI(j)').^2

% replace interface nodes in EToVW with tet interface nodes
EToVTW = EToVW;
for e = 1:num_interface_elems
    v = EToVW(e,:);
    [~, iv, ij] = intersect(v,j);
    v = v + length(VXT) - num_interface_nodes; % increment to include tet nodes
    v(iv) = interface_ids(i(ij)); % swap out wedge interface nodes for tet interface nodes
    EToVTW(e,:) = v;
end
EToVTW(num_interface_elems+1:end,:) = EToVTW(num_interface_elems+1:end,:) + length(VXT) - num_interface_nodes;
EToV = [[EToVT -ones(size(EToVT,1),2)];EToVTW];

% % combine nodes, trim one copy of interface nodes
VX = [VXT(:)' VXW(num_interface_nodes+1:end)];
VY = [VYT(:)' VYW(num_interface_nodes+1:end)];
VZ = [VZT(:)' VZW(num_interface_nodes+1:end)]; 

% drawTetWedgeMesh(VX,VY,VZ,EToV)
% drawTetWedgeMesh(VXW,VYW,VZW,EToVW(1:num_interface_elems,:))
% return
%% write to GMSH file

write_gmsh_file(filename,VX,VY,VZ,EToV)



