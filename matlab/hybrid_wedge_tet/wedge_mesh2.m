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
0; ...
1; ...
1; ...
-1; ...
0; ...
-1; ...
-0.163975; ...
-1; ...
1; ...
0; ...
1; ...
1; ...
-1; ...
0; ...
-1; ...
-0.163975; ...
-1; ...
1; ...
0; ...
1; ...
1; ...
-1; ...
0; ...
-1; ...
-0.163975; ...
];
VY = [ ...
-1; ...
-1; ...
-1; ...
1; ...
0; ...
1; ...
1; ...
0; ...
-0.163975; ...
-1; ...
-1; ...
-1; ...
1; ...
0; ...
1; ...
1; ...
0; ...
-0.163975; ...
-1; ...
-1; ...
-1; ...
1; ...
0; ...
1; ...
1; ...
0; ...
-0.163975; ...
];
VZ = [ ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
0; ...
0; ...
0; ...
0; ...
0; ...
0; ...
0; ...
0; ...
0; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
];
EToV = [ ...
8,...
9,...
7,...
4,...
2,...
7,...
5,...
3,...
17,...
18,...
16,...
13,...
11,...
16,...
14,...
12;...
1,...
1,...
6,...
7,...
5,...
8,...
7,...
5,...
10,...
10,...
15,...
16,...
14,...
17,...
16,...
14;...
9,...
3,...
8,...
5,...
3,...
9,...
9,...
9,...
18,...
12,...
17,...
14,...
12,...
18,...
18,...
18;...
17,...
18,...
16,...
13,...
11,...
16,...
14,...
12,...
26,...
27,...
25,...
22,...
20,...
25,...
23,...
21;...
10,...
10,...
15,...
16,...
14,...
17,...
16,...
14,...
19,...
19,...
24,...
25,...
23,...
26,...
25,...
23;...
18,...
12,...
17,...
14,...
12,...
18,...
18,...
18,...
27,...
21,...
26,...
23,...
21,...
27,...
27,...
27;...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0;...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0;...
];
EToF = [ ...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2;...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0;...
0,...
5,...
0,...
0,...
0,...
4,...
5,...
5,...
0,...
5,...
0,...
0,...
0,...
4,...
5,...
5;...
5,...
4,...
3,...
0,...
0,...
5,...
5,...
4,...
5,...
4,...
3,...
0,...
0,...
5,...
5,...
4;...
3,...
0,...
0,...
3,...
3,...
4,...
4,...
4,...
3,...
0,...
0,...
3,...
3,...
4,...
4,...
4;...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0;...
];
EToE = [ ...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
1,...
2,...
3,...
4,...
5,...
6,...
7,...
8;...
9,...
10,...
11,...
12,...
13,...
14,...
15,...
16,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0;...
0,...
1,...
0,...
0,...
0,...
3,...
4,...
5,...
0,...
9,...
0,...
0,...
0,...
11,...
12,...
13;...
6,...
8,...
6,...
0,...
0,...
7,...
8,...
2,...
14,...
16,...
14,...
0,...
0,...
15,...
16,...
10;...
2,...
0,...
0,...
7,...
8,...
1,...
6,...
7,...
10,...
0,...
0,...
15,...
16,...
9,...
14,...
15;...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0;...
];
EToV = EToV';
EToE = EToE';
EToF = EToF';