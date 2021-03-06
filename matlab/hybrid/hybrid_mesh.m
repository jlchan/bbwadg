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
-1; ...
1; ...
1; ...
-1; ...
1; ...
1; ...
-1; ...
-1; ...
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
-1; ...
-1; ...
1; ...
1; ...
1; ...
-1; ...
-1; ...
1; ...
];
VZ = [ ...
-1; ...
-1; ...
-1; ...
-1; ...
0.33; ...
0.33; ...
0.33; ...
0.33; ...
-0.33; ...
-0.33; ...
-0.33; ...
-0.33; ...
1; ...
1; ...
1; ...
1; ...
];
EToV = [ ...
9,...
10,...
11,...
12,...
7,...
6,...
5,...
8;...
5,...
6,...
13,...
8,...
7,...
16,...
-1,...
-1;...
7,...
16,...
15,...
6,...
13,...
14,...
-1,...
-1;...
1,...
2,...
3,...
4,...
12,...
-1,...
-1,...
-1;...
9,...
12,...
11,...
10,...
2,...
-1,...
-1,...
-1;...
1,...
9,...
2,...
12,...
-1,...
-1,...
-1,...
-1;...
2,...
3,...
12,...
11,...
-1,...
-1,...
-1,...
-1;...
];
EToF = [ ...
5,...
0,...
0,...
0,...
0,...
3;...
0,...
0,...
6,...
0,...
3,...
0;...
0,...
0,...
5,...
0,...
0,...
0;...
3,...
0,...
1,...
0,...
0,...
0;...
4,...
0,...
3,...
0,...
1,...
0;...
0,...
0,...
1,...
1,...
0,...
0;...
3,...
0,...
3,...
0,...
0,...
0;...
];
EToE = [ ...
5,...
0,...
0,...
0,...
0,...
2;...
0,...
0,...
1,...
0,...
3,...
0;...
0,...
0,...
2,...
0,...
0,...
0;...
6,...
0,...
7,...
0,...
0,...
0;...
6,...
0,...
7,...
0,...
1,...
0;...
0,...
0,...
4,...
5,...
0,...
0;...
4,...
0,...
5,...
0,...
0,...
0;...
];

% VZ = VZ*3; % to get square elements

% VXcopy = VX; VX = VZ; VZ = VXcopy;