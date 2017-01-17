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
0; ...
-1; ...
0; ...
-1; ...
1; ...
-1; ...
0; ...
0; ...
-1; ...
0; ...
0; ...
1; ...
0; ...
1; ...
1; ...
0; ...
1; ...
0; ...
-1; ...
-0.5; ...
-1; ...
-0.5; ...
-1; ...
0; ...
-1; ...
-0.5; ...
-0.5; ...
-1; ...
-0.5; ...
-0.5; ...
-1; ...
-1; ...
0; ...
-1; ...
-0.5; ...
-0.5; ...
-1; ...
-0.5; ...
-0.5; ...
-1; ...
-0.5; ...
-0.5; ...
-1; ...
0; ...
-1; ...
-0.5; ...
-0.5; ...
-1; ...
-0.5; ...
-0.5; ...
-1; ...
0; ...
-1; ...
-0.5; ...
-0.5; ...
-1; ...
-0.5; ...
-0.5; ...
0; ...
0.5; ...
0.5; ...
1; ...
0; ...
0.5; ...
0.5; ...
0; ...
0.5; ...
0.5; ...
1; ...
0; ...
0.5; ...
0.5; ...
0; ...
0.5; ...
0.5; ...
1; ...
0.5; ...
1; ...
0.5; ...
1; ...
0; ...
1; ...
0.5; ...
0.5; ...
1; ...
0.5; ...
0.5; ...
1; ...
1; ...
0; ...
1; ...
0.5; ...
0.5; ...
1; ...
0.5; ...
0.5; ...
1; ...
0.5; ...
1; ...
0.5; ...
1; ...
0.5; ...
1; ...
0.5; ...
1; ...
0.5; ...
0; ...
-0.5; ...
0; ...
-0.5; ...
0; ...
-0.5; ...
0; ...
-0.5; ...
0; ...
-0.5; ...
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
-1; ...
-1; ...
0; ...
1; ...
-1; ...
0; ...
-1; ...
0; ...
0; ...
-1; ...
0; ...
0; ...
1; ...
1; ...
0; ...
1; ...
0; ...
1; ...
-0.5; ...
-1; ...
-1; ...
-0.5; ...
0; ...
-1; ...
-0.5; ...
-1; ...
-0.5; ...
-0.5; ...
-1; ...
-0.5; ...
-1; ...
0; ...
-1; ...
-0.5; ...
-1; ...
-0.5; ...
-0.5; ...
-1; ...
-0.5; ...
0.5; ...
0; ...
0.5; ...
1; ...
0; ...
0.5; ...
0; ...
0.5; ...
0.5; ...
0; ...
0.5; ...
1; ...
0; ...
0.5; ...
0; ...
0.5; ...
0.5; ...
0; ...
0.5; ...
-0.5; ...
-1; ...
-0.5; ...
-1; ...
-0.5; ...
-1; ...
-0.5; ...
-0.5; ...
-1; ...
-0.5; ...
-1; ...
-0.5; ...
-1; ...
-0.5; ...
-0.5; ...
-1; ...
-0.5; ...
0.5; ...
1; ...
1; ...
0.5; ...
0; ...
1; ...
0.5; ...
1; ...
0.5; ...
0.5; ...
1; ...
0.5; ...
1; ...
0; ...
1; ...
0.5; ...
1; ...
0.5; ...
0.5; ...
1; ...
0.5; ...
-0.5; ...
0; ...
-0.5; ...
0; ...
-0.5; ...
0; ...
-0.5; ...
0; ...
-0.5; ...
0; ...
0.5; ...
1; ...
0.5; ...
1; ...
0.5; ...
1; ...
0.5; ...
1; ...
0.5; ...
1; ...
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
1; ...
1; ...
0; ...
1; ...
0; ...
0; ...
-1; ...
-1; ...
-1; ...
0; ...
0; ...
0; ...
1; ...
1; ...
0; ...
-1; ...
-1; ...
0; ...
0; ...
1; ...
1; ...
0.5; ...
1; ...
0.5; ...
0.5; ...
0; ...
0; ...
0; ...
0.5; ...
0.5; ...
0.5; ...
-0.5; ...
-0.5; ...
-0.5; ...
-1; ...
-1; ...
-1; ...
-0.5; ...
-0.5; ...
-0.5; ...
1; ...
1; ...
1; ...
0.5; ...
0.5; ...
0; ...
0; ...
0; ...
0.5; ...
0.5; ...
0.5; ...
-0.5; ...
-0.5; ...
-1; ...
-1; ...
-1; ...
-0.5; ...
-0.5; ...
-0.5; ...
1; ...
1; ...
1; ...
0.5; ...
0; ...
0; ...
0; ...
0.5; ...
0.5; ...
0.5; ...
-0.5; ...
-1; ...
-1; ...
-1; ...
-0.5; ...
-0.5; ...
-0.5; ...
1; ...
1; ...
0.5; ...
1; ...
0.5; ...
0.5; ...
0; ...
0; ...
0; ...
0.5; ...
0.5; ...
0.5; ...
-0.5; ...
-0.5; ...
-0.5; ...
-1; ...
-1; ...
-1; ...
-0.5; ...
-0.5; ...
-0.5; ...
1; ...
1; ...
0; ...
0; ...
0.5; ...
0.5; ...
-1; ...
-1; ...
-0.5; ...
-0.5; ...
1; ...
1; ...
0; ...
0; ...
0.5; ...
0.5; ...
-1; ...
-1; ...
-0.5; ...
-0.5; ...
];
EToV = [ ...
7,...
28,...
29,...
30,...
37,...
38,...
-1,...
-1;...
30,...
37,...
38,...
11,...
34,...
35,...
-1,...
-1;...
28,...
9,...
31,...
37,...
32,...
39,...
-1,...
-1;...
37,...
32,...
39,...
34,...
18,...
36,...
-1,...
-1;...
29,...
31,...
10,...
38,...
39,...
33,...
-1,...
-1;...
38,...
39,...
33,...
35,...
36,...
19,...
-1,...
-1;...
31,...
29,...
28,...
39,...
38,...
37,...
-1,...
-1;...
39,...
38,...
37,...
36,...
35,...
34,...
-1,...
-1;...
11,...
34,...
35,...
40,...
46,...
47,...
-1,...
-1;...
40,...
46,...
47,...
1,...
43,...
44,...
-1,...
-1;...
34,...
18,...
36,...
46,...
41,...
48,...
-1,...
-1;...
46,...
41,...
48,...
43,...
15,...
45,...
-1,...
-1;...
35,...
36,...
19,...
47,...
48,...
42,...
-1,...
-1;...
47,...
48,...
42,...
44,...
45,...
16,...
-1,...
-1;...
36,...
35,...
34,...
48,...
47,...
46,...
-1,...
-1;...
48,...
47,...
46,...
45,...
44,...
43,...
-1,...
-1;...
9,...
49,...
50,...
32,...
57,...
58,...
-1,...
-1;...
32,...
57,...
58,...
18,...
54,...
55,...
-1,...
-1;...
49,...
8,...
51,...
57,...
52,...
59,...
-1,...
-1;...
57,...
52,...
59,...
54,...
13,...
56,...
-1,...
-1;...
50,...
51,...
12,...
58,...
59,...
53,...
-1,...
-1;...
58,...
59,...
53,...
55,...
56,...
20,...
-1,...
-1;...
51,...
50,...
49,...
59,...
58,...
57,...
-1,...
-1;...
59,...
58,...
57,...
56,...
55,...
54,...
-1,...
-1;...
18,...
54,...
55,...
41,...
65,...
66,...
-1,...
-1;...
41,...
65,...
66,...
15,...
62,...
63,...
-1,...
-1;...
54,...
13,...
56,...
65,...
60,...
67,...
-1,...
-1;...
65,...
60,...
67,...
62,...
4,...
64,...
-1,...
-1;...
55,...
56,...
20,...
66,...
67,...
61,...
-1,...
-1;...
66,...
67,...
61,...
63,...
64,...
17,...
-1,...
-1;...
56,...
55,...
54,...
67,...
66,...
65,...
-1,...
-1;...
67,...
66,...
65,...
64,...
63,...
62,...
-1,...
-1;...
10,...
68,...
69,...
33,...
75,...
76,...
-1,...
-1;...
33,...
75,...
76,...
19,...
72,...
73,...
-1,...
-1;...
68,...
12,...
70,...
75,...
53,...
77,...
-1,...
-1;...
75,...
53,...
77,...
72,...
20,...
74,...
-1,...
-1;...
69,...
70,...
6,...
76,...
77,...
71,...
-1,...
-1;...
76,...
77,...
71,...
73,...
74,...
14,...
-1,...
-1;...
70,...
69,...
68,...
77,...
76,...
75,...
-1,...
-1;...
77,...
76,...
75,...
74,...
73,...
72,...
-1,...
-1;...
19,...
72,...
73,...
42,...
82,...
83,...
-1,...
-1;...
42,...
82,...
83,...
16,...
79,...
80,...
-1,...
-1;...
72,...
20,...
74,...
82,...
61,...
84,...
-1,...
-1;...
82,...
61,...
84,...
79,...
17,...
81,...
-1,...
-1;...
73,...
74,...
14,...
83,...
84,...
78,...
-1,...
-1;...
83,...
84,...
78,...
80,...
81,...
2,...
-1,...
-1;...
74,...
73,...
72,...
84,...
83,...
82,...
-1,...
-1;...
84,...
83,...
82,...
81,...
80,...
79,...
-1,...
-1;...
12,...
68,...
50,...
53,...
75,...
58,...
-1,...
-1;...
53,...
75,...
58,...
20,...
72,...
55,...
-1,...
-1;...
68,...
10,...
31,...
75,...
33,...
39,...
-1,...
-1;...
75,...
33,...
39,...
72,...
19,...
36,...
-1,...
-1;...
50,...
31,...
9,...
58,...
39,...
32,...
-1,...
-1;...
58,...
39,...
32,...
55,...
36,...
18,...
-1,...
-1;...
31,...
50,...
68,...
39,...
58,...
75,...
-1,...
-1;...
39,...
58,...
75,...
36,...
55,...
72,...
-1,...
-1;...
20,...
72,...
55,...
61,...
82,...
66,...
-1,...
-1;...
61,...
82,...
66,...
17,...
79,...
63,...
-1,...
-1;...
72,...
19,...
36,...
82,...
42,...
48,...
-1,...
-1;...
82,...
42,...
48,...
79,...
16,...
45,...
-1,...
-1;...
55,...
36,...
18,...
66,...
48,...
41,...
-1,...
-1;...
66,...
48,...
41,...
63,...
45,...
15,...
-1,...
-1;...
36,...
55,...
72,...
48,...
66,...
82,...
-1,...
-1;...
48,...
66,...
82,...
45,...
63,...
79,...
-1,...
-1;...
5,...
85,...
86,...
87,...
94,...
95,...
-1,...
-1;...
87,...
94,...
95,...
23,...
91,...
92,...
-1,...
-1;...
85,...
21,...
88,...
94,...
89,...
96,...
-1,...
-1;...
94,...
89,...
96,...
91,...
26,...
93,...
-1,...
-1;...
86,...
88,...
22,...
95,...
96,...
90,...
-1,...
-1;...
95,...
96,...
90,...
92,...
93,...
27,...
-1,...
-1;...
88,...
86,...
85,...
96,...
95,...
94,...
-1,...
-1;...
96,...
95,...
94,...
93,...
92,...
91,...
-1,...
-1;...
23,...
91,...
92,...
97,...
103,...
104,...
-1,...
-1;...
97,...
103,...
104,...
3,...
100,...
101,...
-1,...
-1;...
91,...
26,...
93,...
103,...
98,...
105,...
-1,...
-1;...
103,...
98,...
105,...
100,...
24,...
102,...
-1,...
-1;...
92,...
93,...
27,...
104,...
105,...
99,...
-1,...
-1;...
104,...
105,...
99,...
101,...
102,...
25,...
-1,...
-1;...
93,...
92,...
91,...
105,...
104,...
103,...
-1,...
-1;...
105,...
104,...
103,...
102,...
101,...
100,...
-1,...
-1;...
21,...
106,...
107,...
89,...
110,...
111,...
-1,...
-1;...
89,...
110,...
111,...
26,...
108,...
109,...
-1,...
-1;...
106,...
6,...
70,...
110,...
71,...
77,...
-1,...
-1;...
110,...
71,...
77,...
108,...
14,...
74,...
-1,...
-1;...
107,...
70,...
12,...
111,...
77,...
53,...
-1,...
-1;...
111,...
77,...
53,...
109,...
74,...
20,...
-1,...
-1;...
70,...
107,...
106,...
77,...
111,...
110,...
-1,...
-1;...
77,...
111,...
110,...
74,...
109,...
108,...
-1,...
-1;...
26,...
108,...
109,...
98,...
114,...
115,...
-1,...
-1;...
98,...
114,...
115,...
24,...
112,...
113,...
-1,...
-1;...
108,...
14,...
74,...
114,...
78,...
84,...
-1,...
-1;...
114,...
78,...
84,...
112,...
2,...
81,...
-1,...
-1;...
109,...
74,...
20,...
115,...
84,...
61,...
-1,...
-1;...
115,...
84,...
61,...
113,...
81,...
17,...
-1,...
-1;...
74,...
109,...
108,...
84,...
115,...
114,...
-1,...
-1;...
84,...
115,...
114,...
81,...
113,...
112,...
-1,...
-1;...
22,...
116,...
117,...
90,...
120,...
121,...
-1,...
-1;...
90,...
120,...
121,...
27,...
118,...
119,...
-1,...
-1;...
116,...
12,...
51,...
120,...
53,...
59,...
-1,...
-1;...
120,...
53,...
59,...
118,...
20,...
56,...
-1,...
-1;...
117,...
51,...
8,...
121,...
59,...
52,...
-1,...
-1;...
121,...
59,...
52,...
119,...
56,...
13,...
-1,...
-1;...
51,...
117,...
116,...
59,...
121,...
120,...
-1,...
-1;...
59,...
121,...
120,...
56,...
119,...
118,...
-1,...
-1;...
27,...
118,...
119,...
99,...
124,...
125,...
-1,...
-1;...
99,...
124,...
125,...
25,...
122,...
123,...
-1,...
-1;...
118,...
20,...
56,...
124,...
61,...
67,...
-1,...
-1;...
124,...
61,...
67,...
122,...
17,...
64,...
-1,...
-1;...
119,...
56,...
13,...
125,...
67,...
60,...
-1,...
-1;...
125,...
67,...
60,...
123,...
64,...
4,...
-1,...
-1;...
56,...
119,...
118,...
67,...
125,...
124,...
-1,...
-1;...
67,...
125,...
124,...
64,...
123,...
122,...
-1,...
-1;...
12,...
116,...
107,...
53,...
120,...
111,...
-1,...
-1;...
53,...
120,...
111,...
20,...
118,...
109,...
-1,...
-1;...
116,...
22,...
88,...
120,...
90,...
96,...
-1,...
-1;...
120,...
90,...
96,...
118,...
27,...
93,...
-1,...
-1;...
107,...
88,...
21,...
111,...
96,...
89,...
-1,...
-1;...
111,...
96,...
89,...
109,...
93,...
26,...
-1,...
-1;...
88,...
107,...
116,...
96,...
111,...
120,...
-1,...
-1;...
96,...
111,...
120,...
93,...
109,...
118,...
-1,...
-1;...
20,...
118,...
109,...
61,...
124,...
115,...
-1,...
-1;...
61,...
124,...
115,...
17,...
122,...
113,...
-1,...
-1;...
118,...
27,...
93,...
124,...
99,...
105,...
-1,...
-1;...
124,...
99,...
105,...
122,...
25,...
102,...
-1,...
-1;...
109,...
93,...
26,...
115,...
105,...
98,...
-1,...
-1;...
115,...
105,...
98,...
113,...
102,...
24,...
-1,...
-1;...
93,...
109,...
118,...
105,...
115,...
124,...
-1,...
-1;...
105,...
115,...
124,...
102,...
113,...
122,...
-1,...
-1;...
];
EToF = [ ...
];
EToE = [ ...
];