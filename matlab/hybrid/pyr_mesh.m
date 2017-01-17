%------------------------------------------------------
% output of mesh3d
% description: c++ codes will read mesh files and
%              write information to matlab script
% output:
% VX,VY,VZ: coordinates of vertices
% EToV: elements to vertices table
% EToE: elements to elements table
% EToV: elements to faces table
%------------------------------------------------------
VX = [ ...
-1; ...
1; ...
1; ...
-1; ...
1; ...
1; ...
-1; ...
-1; ...
0; ...
];
VY = [ ...
-1; ...
-1; ...
1; ...
1; ...
1; ...
-1; ...
-1; ...
1; ...
0; ...
];
VZ = [ ...
-1; ...
-1; ...
-1; ...
-1; ...
1; ...
1; ...
1; ...
1; ...
0; ...
];
EToV = [ ...
1,...
2,...
3,...
4,...
9,...
-1,...
-1,...
-1;...
4,...
3,...
5,...
8,...
9,...
-1,...
-1,...
-1;...
5,...
6,...
7,...
8,...
9,...
-1,...
-1,...
-1;...
2,...
1,...
7,...
6,...
9,...
-1,...
-1,...
-1;...
2,...
6,...
5,...
3,...
9,...
-1,...
-1,...
-1;...
1,...
4,...
8,...
7,...
9,...
-1,...
-1,...
-1;...
];
EToF = [ ...
1,...
1,...
2,...
1,...
0,...
0;...
4,...
3,...
4,...
2,...
0,...
0;...
3,...
4,...
4,...
4,...
0,...
0;...
1,...
1,...
2,...
3,...
0,...
0;...
2,...
3,...
1,...
3,...
0,...
0;...
2,...
3,...
2,...
4,...
0,...
0;...
];
EToE = [ ...
4,...
6,...
5,...
2,...
0,...
0;...
1,...
6,...
5,...
3,...
0,...
0;...
5,...
2,...
4,...
6,...
0,...
0;...
1,...
5,...
6,...
3,...
0,...
0;...
4,...
1,...
3,...
2,...
0,...
0;...
1,...
4,...
2,...
3,...
0,...
0;...
];
