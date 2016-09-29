
Globals2D;

N = 1;

%filename = 'Grid/Other/circA01.neu';
% filename = 'Grid/CFD/cylinderA00075c.neu';
% filename = 'Grid/CFD/cylinderCA0015.neu';
% filename = 'Grid/CNS2D/CNScylK930.neu';

% filename = 'Grid/CNS2D/bumpSA025.neu';
filename = 'Grid/CNS2D/couette.neu';
% filename = 'Grid/CNS2D/cylSA02.neu';

% filename = 'Grid/Euler2D/inlet1K546.neu';
% filename = 'Grid/CFD/pvortex3A025.neu';
[Nv, VX, VY, K, EToV, BCType] = MeshReaderGambitBC2D(filename);

% This builds the nodal DG stuff
StartUp2D;
BuildBCMaps2D;

PlotMesh2D