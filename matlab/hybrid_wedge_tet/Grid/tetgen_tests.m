%"/Users/jlchan/Desktop/total/matlab/iso2mesh/bin/tetgen.mexmaci64" -A -q1.414a25  "/private/tmp/iso2mesh-jlchan/post_vmesh.poly"

clear
maxvol = .25;
file = '"Grid/tetgen/cube.poly"';
file = '"Grid/tetgen/layer_mesh.poly"';
system(['"/Users/jlchan/Desktop/total/matlab/iso2mesh/bin/tetgen.mexmaci64" -pq1.414Va' num2str(maxvol) ' ' file])

% read in the generated mesh
% [node,elem,face]=readtetgen('Grid/tetgen/cube.1');
[node,elem,face]=readtetgen('Grid/tetgen/layer_mesh.1');
% trisurf(face(:,1:3),node(:,1),node(:,2),node(:,3));
plotmesh(node(:,1:3),elem);
% axis equal;

% fprintf(1,'volume mesh generation is complete\n');

