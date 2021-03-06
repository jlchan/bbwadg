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
-1; ...
1; ...
0; ...
0; ...
-0.5; ...
-1; ...
0; ...
-0.416667; ...
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
0; ...
0; ...
-1; ...
1; ...
-1; ...
-0.5; ...
0; ...
-0.416667; ...
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
0.057462; ...
-0.117163; ...
0.071731; ...
-0.052572; ...
0.085725; ...
0.025784; ...
0.032186; ...
-0.042044; ...
0.027444; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
1.02124; ...
1.02124; ...
1.02124; ...
1.02124; ...
1.02124; ...
1.02124; ...
1.02124; ...
1.02124; ...
1.02124; ...
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
19,...
18,...
7,...
5,...
15,...
12,...
6,...
5,...
21,...
17,...
9,...
13,...
12,...
18,...
19,...
16,...
7,...
13,...
17,...
8,...
14,...
12,...
16,...
11,...
15,...
6,...
13,...
11,...
9,...
11,...
9;...
5,...
5,...
10,...
11,...
9,...
12,...
11,...
9,...
14,...
5,...
16,...
18,...
2,...
21,...
2,...
12,...
12,...
12,...
6,...
16,...
10,...
1,...
12,...
6,...
5,...
20,...
3,...
3,...
4,...
4,...
13,...
17,...
6,...
13,...
12,...
12,...
13,...
8,...
8;...
13,...
7,...
12,...
9,...
7,...
13,...
13,...
13,...
12,...
7,...
18,...
19,...
6,...
20,...
16,...
13,...
20,...
11,...
7,...
7,...
17,...
5,...
5,...
13,...
13,...
16,...
8,...
15,...
12,...
17,...
6,...
8,...
9,...
9,...
11,...
17,...
11,...
9,...
15;...
29,...
30,...
28,...
25,...
23,...
28,...
26,...
24,...
21,...
21,...
21,...
21,...
20,...
13,...
20,...
21,...
14,...
10,...
13,...
21,...
4,...
19,...
21,...
7,...
21,...
21,...
20,...
20,...
20,...
20,...
20,...
20,...
20,...
20,...
20,...
20,...
20,...
20,...
20;...
22,...
22,...
27,...
28,...
26,...
29,...
28,...
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
30,...
24,...
29,...
26,...
24,...
30,...
30,...
30,...
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
1,...
1,...
4,...
1,...
1,...
1,...
1,...
3,...
0,...
0,...
0,...
3,...
0,...
1,...
0,...
1,...
1,...
1,...
1,...
3,...
2,...
0,...
0,...
1,...
1,...
2,...
0,...
0,...
0,...
3,...
1,...
0,...
0,...
2,...
1,...
1,...
1,...
1,...
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
2,...
4,...
2,...
0,...
4,...
4,...
4,...
4,...
1,...
1,...
3,...
0,...
0,...
3,...
0,...
4,...
4,...
0,...
4,...
0,...
4,...
1,...
3,...
3,...
4,...
3,...
4,...
4,...
3,...
4;...
0,...
5,...
0,...
0,...
0,...
4,...
5,...
5,...
2,...
3,...
3,...
3,...
2,...
2,...
3,...
4,...
0,...
0,...
1,...
3,...
1,...
1,...
3,...
1,...
3,...
2,...
4,...
4,...
4,...
4,...
3,...
2,...
3,...
4,...
4,...
2,...
3,...
3,...
3;...
5,...
4,...
3,...
0,...
0,...
5,...
5,...
4,...
2,...
2,...
0,...
0,...
2,...
2,...
0,...
2,...
3,...
1,...
4,...
2,...
0,...
0,...
2,...
4,...
3,...
0,...
2,...
0,...
2,...
0,...
2,...
3,...
3,...
2,...
2,...
3,...
3,...
2,...
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
16,...
25,...
18,...
38,...
19,...
35,...
37,...
19,...
0,...
0,...
0,...
22,...
0,...
17,...
0,...
1,...
14,...
36,...
5,...
24,...
18,...
0,...
0,...
31,...
2,...
31,...
0,...
0,...
0,...
21,...
24,...
0,...
0,...
19,...
6,...
18,...
7,...
4,...
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
12,...
20,...
10,...
0,...
16,...
13,...
23,...
9,...
21,...
34,...
26,...
0,...
0,...
9,...
0,...
10,...
14,...
0,...
27,...
0,...
29,...
26,...
36,...
13,...
31,...
14,...
35,...
34,...
32,...
38;...
0,...
1,...
0,...
0,...
0,...
3,...
4,...
5,...
23,...
11,...
10,...
23,...
33,...
35,...
31,...
25,...
0,...
0,...
8,...
25,...
30,...
12,...
12,...
20,...
20,...
20,...
32,...
39,...
17,...
36,...
15,...
38,...
39,...
33,...
37,...
32,...
38,...
37,...
33;...
6,...
8,...
6,...
0,...
0,...
7,...
8,...
2,...
17,...
25,...
0,...
0,...
15,...
26,...
0,...
14,...
29,...
3,...
24,...
11,...
0,...
0,...
16,...
19,...
16,...
0,...
28,...
0,...
30,...
0,...
34,...
27,...
34,...
37,...
36,...
30,...
35,...
39,...
28;...
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
