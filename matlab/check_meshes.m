
Globals2D;

N = 1;

filename = 'Grid/Other/circleK8.neu';
% filename = 'Grid/CFD/cylinderA00075c.neu';
% filename = 'Grid/CFD/cylinderCA0015.neu';
% filename = 'Grid/CNS2D/CNScylK930.neu';

% filename = 'Grid/CNS2D/bumpSA025.neu';
filename = 'Grid/CNS2D/couette.neu';
% filename = 'Grid/CNS2D/cylSA02.neu';

% filename = 'Grid/Euler2D/inlet1K546.neu';
% filename = 'Grid/CFD/pvortex3A025.neu';
% filename = 'Grid/CFD/kovA02.neu';
% filename = 'Grid/Euler2D/forward3A00004.neu';
[Nv, VY, VX, K, EToV, BCType] = MeshReaderGambitBC2D(filename);
VX = VX - .5;
VY = VY - .25;
% This builds the nodal DG stuff
StartUp2D;
BuildBCMaps2D;

PlotMesh2D
axis on