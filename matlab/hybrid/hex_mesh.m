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

% VX = [-1 0 1 -1 0 1 -1 0 1 ...
%     -1 0 1 -1 0 1 -1 0 1 ...
%     -1 0 1 -1 0 1 -1 0 1];
% 
% VY = [-1 -1 -1 0 0 0 1 1 1 ...
%     -1 -1 -1 0 0 0 1 1 1 ...
%     -1 -1 -1 0 0 0 1 1 1];
% 
% VZ = [-ones(1,9) zeros(1,9) ones(1,9)];
    
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
0; ...
-1; ...
-1; ...
1; ...
1; ...
0; ...
1; ...
-1; ...
0; ...
-1; ...
1; ...
0; ...
0; ...
0; ...
-1; ...
1; ...
0; ...
0; ...
0; ...
0; ...
-1; ...
-1; ...
1; ...
1; ...
0; ...
1; ...
-1; ...
0; ...
0; ...
-1; ...
1; ...
0; ...
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
-1; ...
-1; ...
1; ...
1; ...
-1; ...
0; ...
-1; ...
0; ...
-1; ...
1; ...
1; ...
1; ...
-1; ...
0; ...
0; ...
1; ...
0; ...
-1; ...
0; ...
0; ...
1; ...
0; ...
0; ...
-1; ...
0; ...
-1; ...
0; ...
-1; ...
1; ...
1; ...
1; ...
0; ...
-1; ...
0; ...
0; ...
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
0; ...
0; ...
0; ...
-1; ...
-1; ...
-0.5; ...
-1; ...
-0.5; ...
-1; ...
-0.5; ...
-0.5; ...
0; ...
0; ...
0; ...
0; ...
-1; ...
-0.5; ...
-0.5; ...
-0.5; ...
-0.5; ...
0; ...
-0.5; ...
1; ...
1; ...
0.5; ...
1; ...
0.5; ...
1; ...
0.5; ...
0.5; ...
1; ...
0.5; ...
0.5; ...
0.5; ...
0.5; ...
0.5; ...
];
EToV = [ ...
1,...
13,...
25,...
14,...
15,...
26,...
31,...
27;...
15,...
26,...
31,...
27,...
9,...
21,...
30,...
22;...
13,...
2,...
16,...
25,...
26,...
17,...
28,...
31;...
26,...
17,...
28,...
31,...
21,...
10,...
23,...
30;...
14,...
25,...
18,...
4,...
27,...
31,...
29,...
20;...
27,...
31,...
29,...
20,...
22,...
30,...
24,...
12;...
25,...
16,...
3,...
18,...
31,...
28,...
19,...
29;...
31,...
28,...
19,...
29,...
30,...
23,...
11,...
24;...
40,...
32,...
7,...
33,...
45,...
41,...
34,...
42;...
45,...
41,...
34,...
42,...
30,...
21,...
9,...
22;...
35,...
6,...
32,...
40,...
43,...
36,...
41,...
45;...
43,...
36,...
41,...
45,...
23,...
10,...
21,...
30;...
37,...
40,...
33,...
8,...
44,...
45,...
42,...
39;...
44,...
45,...
42,...
39,...
24,...
30,...
22,...
12;...
5,...
35,...
40,...
37,...
38,...
43,...
45,...
44;...
38,...
43,...
45,...
44,...
11,...
23,...
30,...
24;...
];
EToF = [ ...
0,...
0,...
5,...
2,...
0,...
1;...
6,...
0,...
5,...
2,...
0,...
6;...
0,...
0,...
0,...
2,...
3,...
1;...
6,...
0,...
0,...
2,...
3,...
6;...
0,...
4,...
5,...
0,...
0,...
1;...
6,...
4,...
5,...
0,...
0,...
6;...
0,...
4,...
0,...
0,...
3,...
1;...
6,...
4,...
0,...
0,...
3,...
6;...
0,...
4,...
0,...
0,...
3,...
1;...
6,...
4,...
0,...
0,...
3,...
6;...
0,...
0,...
0,...
2,...
3,...
1;...
6,...
0,...
0,...
2,...
3,...
6;...
0,...
4,...
5,...
0,...
0,...
1;...
6,...
4,...
5,...
0,...
0,...
6;...
0,...
0,...
5,...
2,...
0,...
1;...
6,...
0,...
5,...
2,...
0,...
6;...
];
EToE = [ ...
0,...
0,...
3,...
5,...
0,...
2;...
1,...
0,...
4,...
6,...
0,...
10;...
0,...
0,...
0,...
7,...
1,...
4;...
3,...
0,...
0,...
8,...
2,...
12;...
0,...
1,...
7,...
0,...
0,...
6;...
5,...
2,...
8,...
0,...
0,...
14;...
0,...
3,...
0,...
0,...
5,...
8;...
7,...
4,...
0,...
0,...
6,...
16;...
0,...
11,...
0,...
0,...
13,...
10;...
9,...
12,...
0,...
0,...
14,...
2;...
0,...
0,...
0,...
9,...
15,...
12;...
11,...
0,...
0,...
10,...
16,...
4;...
0,...
15,...
9,...
0,...
0,...
14;...
13,...
16,...
10,...
0,...
0,...
6;...
0,...
0,...
11,...
13,...
0,...
16;...
15,...
0,...
12,...
14,...
0,...
8;...
];