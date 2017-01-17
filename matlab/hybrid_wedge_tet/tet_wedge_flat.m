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
-1; ...
1; ...
0; ...
1; ...
1; ...
-1; ...
0; ...
-1; ...
-0.163975; ...
0; ...
1; ...
0; ...
-1; ...
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
1; ...
1; ...
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
0; ...
1; ...
0; ...
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
0; ...
0; ...
0; ...
0; ...
0; ...
0; ...
0; ...
0; ...
0; ...
-1; ...
-1; ...
-1; ...
-1; ...
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
12,...
13,...
11,...
8,...
6,...
11,...
9,...
7,...
9,...
13,...
16,...
7,...
16,...
13,...
7,...
9,...
9,...
4,...
16,...
16,...
13,...
17,...
17,...
17,...
15,...
13,...
16,...
11,...
14,...
1;...
5,...
5,...
10,...
11,...
9,...
12,...
11,...
9,...
6,...
1,...
11,...
9,...
9,...
12,...
14,...
7,...
11,...
12,...
9,...
14,...
14,...
11,...
13,...
11,...
7,...
15,...
9,...
16,...
13,...
12;...
13,...
7,...
12,...
9,...
7,...
13,...
13,...
13,...
7,...
7,...
17,...
15,...
3,...
1,...
1,...
15,...
8,...
11,...
11,...
15,...
15,...
12,...
14,...
4,...
14,...
16,...
15,...
17,...
17,...
13;...
25,...
26,...
24,...
21,...
19,...
24,...
22,...
20,...
2,...
5,...
13,...
13,...
11,...
5,...
13,...
2,...
3,...
10,...
13,...
13,...
7,...
13,...
16,...
12,...
2,...
9,...
3,...
4,...
1,...
17;...
18,...
18,...
23,...
24,...
22,...
25,...
24,...
22,...
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
0,...
0,...
0,...
0,...
0,...
0,...
0;...
26,...
20,...
25,...
22,...
20,...
26,...
26,...
26,...
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
2,...
3,...
4,...
1,...
1,...
4,...
4,...
2,...
1,...
3,...
1,...
1,...
2,...
1,...
0,...
1,...
1,...
4,...
2,...
0,...
4,...
2,...
1,...
3,...
4,...
3,...
4,...
1,...
1,...
1;...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
3,...
3,...
1,...
1,...
1,...
2,...
3,...
4,...
0,...
3,...
4,...
2,...
4,...
3,...
1,...
4,...
4,...
1,...
0,...
4,...
0;...
0,...
5,...
0,...
0,...
0,...
4,...
5,...
5,...
2,...
1,...
2,...
3,...
0,...
2,...
1,...
0,...
0,...
0,...
2,...
1,...
3,...
4,...
0,...
0,...
0,...
2,...
0,...
1,...
0,...
4;...
5,...
4,...
3,...
0,...
0,...
5,...
5,...
4,...
0,...
0,...
2,...
2,...
2,...
0,...
2,...
2,...
0,...
1,...
1,...
1,...
1,...
1,...
2,...
1,...
0,...
1,...
0,...
0,...
3,...
3;...
3,...
0,...
0,...
3,...
3,...
4,...
4,...
4,...
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
14,...
10,...
18,...
17,...
9,...
22,...
19,...
12,...
5,...
15,...
28,...
16,...
27,...
30,...
0,...
12,...
4,...
24,...
13,...
0,...
20,...
24,...
29,...
28,...
21,...
20,...
26,...
11,...
23,...
14;...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
14,...
19,...
8,...
19,...
1,...
21,...
9,...
13,...
0,...
26,...
23,...
15,...
11,...
11,...
22,...
16,...
12,...
13,...
0,...
15,...
0;...
0,...
1,...
0,...
0,...
0,...
3,...
4,...
5,...
16,...
2,...
23,...
21,...
0,...
10,...
10,...
0,...
0,...
0,...
11,...
26,...
12,...
30,...
0,...
0,...
0,...
19,...
0,...
24,...
0,...
29;...
6,...
8,...
6,...
0,...
0,...
7,...
8,...
2,...
0,...
0,...
22,...
26,...
17,...
0,...
29,...
25,...
0,...
3,...
7,...
21,...
25,...
6,...
20,...
18,...
0,...
27,...
0,...
0,...
30,...
22;...
2,...
0,...
0,...
7,...
8,...
1,...
6,...
7,...
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
