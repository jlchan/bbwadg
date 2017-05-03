% function inclusion

%  Here are all the nodes to be used.
%
node = [
    0.0, 0.0;
    3.0, 0.0;
    3.0, 3.0;
    0.0, 3.0;
    1.0, 1.5;
    2.0, 1.0;
    2.0, 2.0];
%
%  Here are all the edges.
%
edge = [
    1,2;
    2,3;
    3,4;
    4,1;
    5,6;
    6,7;
    7,5];

%  Face 1 is the outer strip.
%  Face 2 is the inner square.

face{1} = [1,2,3,4,5,6,7];
face{2} = [5,6,7];

%  HDATA allows us to set options.
%    HDATA.HMAX is the maximum triangle side length.

hdata = [];
hdata.hmax = .1;
%
%  MESHFACES output is:
%
%    P     = Nx2 array of nodal XY co-ordinates.
%    T     = Mx3 array of triangles as indicies into P, defined with a
%            counter-clockwise node ordering.
%    FNUM  = Mx1 array of face numbers for each triangle in T.
%
[ p, t, fnum, stats ] = meshfaces ( node, edge, face, hdata );

axis on
