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
0.010732; ...
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
-0.482763; ...
-1; ...
-0.5; ...
-0.48636; ...
-1; ...
-1; ...
0; ...
-1; ...
-0.5; ...
-0.5; ...
-1; ...
-0.5; ...
-0.516868; ...
-1; ...
-0.5; ...
-0.5; ...
-1; ...
0.010268; ...
-1; ...
-0.484317; ...
-0.488716; ...
-1; ...
-0.510346; ...
-0.487622; ...
-1; ...
-0.012862; ...
-1; ...
-0.5; ...
-0.5; ...
-1; ...
-0.513673; ...
-0.48434; ...
0; ...
0.5; ...
0.5; ...
1; ...
0.01122; ...
0.5; ...
0.512553; ...
0.011913; ...
0.5; ...
0.517195; ...
1; ...
0; ...
0.5; ...
0.5; ...
-0.011931; ...
0.5; ...
0.511826; ...
1; ...
0.5; ...
1; ...
0.5; ...
1; ...
0; ...
1; ...
0.5; ...
0.488117; ...
1; ...
0.5; ...
0.483769; ...
1; ...
1; ...
0; ...
1; ...
0.5; ...
0.5; ...
1; ...
0.5; ...
0.512144; ...
1; ...
0.5; ...
1; ...
0.485741; ...
1; ...
0.515844; ...
1; ...
0.5; ...
1; ...
0.51352; ...
0; ...
-0.5; ...
-0.015957; ...
-0.5; ...
-0.014515; ...
-0.5; ...
0; ...
-0.5; ...
0.015611; ...
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
0.012089; ...
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
-0.488818; ...
-0.5; ...
-1; ...
-0.483998; ...
-1; ...
0; ...
-1; ...
-0.5; ...
-1; ...
-0.5; ...
-0.5; ...
-1; ...
-0.515942; ...
0.5; ...
0; ...
0.5; ...
1; ...
0.016368; ...
0.5; ...
0.015573; ...
0.515295; ...
0.5; ...
-0.010728; ...
0.517127; ...
1; ...
-0.015741; ...
0.5; ...
0; ...
0.5; ...
0.5; ...
-0.013342; ...
0.51207; ...
-0.5; ...
-1; ...
-0.5; ...
-1; ...
-0.489108; ...
-1; ...
-0.48561; ...
-0.486205; ...
-1; ...
-0.485896; ...
-1; ...
-0.5; ...
-1; ...
-0.5; ...
-0.516305; ...
-1; ...
-0.483031; ...
0.5; ...
1; ...
1; ...
0.5; ...
0; ...
1; ...
0.5; ...
1; ...
0.48538; ...
0.5; ...
1; ...
0.485611; ...
1; ...
0; ...
1; ...
0.5; ...
1; ...
0.5; ...
0.5; ...
1; ...
0.515679; ...
-0.5; ...
0; ...
-0.5; ...
-0.010569; ...
-0.5; ...
0.017005; ...
-0.5; ...
0; ...
-0.5; ...
0.010089; ...
0.5; ...
1; ...
0.487666; ...
1; ...
0.488028; ...
1; ...
0.5; ...
1; ...
0.513379; ...
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
0.014102; ...
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
0.017279; ...
0.5; ...
0.5; ...
0.511064; ...
-0.5; ...
-0.5; ...
-0.5; ...
-1; ...
-1; ...
-1; ...
-0.5; ...
-0.5; ...
-0.517196; ...
1; ...
1; ...
1; ...
0.5; ...
0.517005; ...
0; ...
0.012942; ...
0.010239; ...
0.5; ...
0.483824; ...
0.510258; ...
-0.5; ...
-0.515964; ...
-1; ...
-1; ...
-1; ...
-0.5; ...
-0.514847; ...
-0.484902; ...
1; ...
1; ...
1; ...
0.5; ...
0.013738; ...
0; ...
0.011679; ...
0.515243; ...
0.5; ...
0.51104; ...
-0.5; ...
-1; ...
-1; ...
-1; ...
-0.511907; ...
-0.5; ...
-0.487375; ...
1; ...
1; ...
0.5; ...
1; ...
0.5; ...
0.5; ...
0; ...
0; ...
-0.01355; ...
0.5; ...
0.5; ...
0.485877; ...
-0.5; ...
-0.5; ...
-0.5; ...
-1; ...
-1; ...
-1; ...
-0.5; ...
-0.5; ...
-0.484347; ...
1; ...
1; ...
0; ...
-0.010405; ...
0.5; ...
0.510974; ...
-1; ...
-1; ...
-0.5; ...
-0.487472; ...
1; ...
1; ...
-0.013964; ...
0; ...
0.485094; ...
0.5; ...
-1; ...
-1; ...
-0.489371; ...
-0.5; ...
];
EToV = [ ...
7,...
30,...
28,...
37,...
29,...
38,...
31,...
39,...
11,...
40,...
34,...
46,...
35,...
47,...
36,...
48,...
9,...
32,...
49,...
57,...
50,...
58,...
51,...
59,...
18,...
41,...
54,...
65,...
55,...
66,...
56,...
67,...
10,...
33,...
68,...
75,...
69,...
76,...
70,...
77,...
19,...
42,...
72,...
82,...
73,...
83,...
74,...
84,...
12,...
53,...
68,...
75,...
50,...
58,...
31,...
39,...
20,...
61,...
72,...
82,...
55,...
66,...
36,...
48,...
5,...
87,...
85,...
94,...
86,...
95,...
88,...
96,...
23,...
97,...
91,...
103,...
92,...
104,...
93,...
105,...
21,...
89,...
106,...
110,...
107,...
111,...
70,...
77,...
26,...
98,...
108,...
114,...
109,...
115,...
74,...
84,...
22,...
90,...
116,...
120,...
117,...
121,...
51,...
59,...
27,...
99,...
118,...
124,...
119,...
125,...
56,...
67,...
12,...
53,...
116,...
120,...
107,...
111,...
88,...
96,...
20,...
61,...
118,...
124,...
109,...
115,...
93,...
105;...
28,...
37,...
9,...
32,...
31,...
39,...
29,...
38,...
34,...
46,...
18,...
41,...
36,...
48,...
35,...
47,...
49,...
57,...
8,...
52,...
51,...
59,...
50,...
58,...
54,...
65,...
13,...
60,...
56,...
67,...
55,...
66,...
68,...
75,...
12,...
53,...
70,...
77,...
69,...
76,...
72,...
82,...
20,...
61,...
74,...
84,...
73,...
83,...
68,...
75,...
10,...
33,...
31,...
39,...
50,...
58,...
72,...
82,...
19,...
42,...
36,...
48,...
55,...
66,...
85,...
94,...
21,...
89,...
88,...
96,...
86,...
95,...
91,...
103,...
26,...
98,...
93,...
105,...
92,...
104,...
106,...
110,...
6,...
71,...
70,...
77,...
107,...
111,...
108,...
114,...
14,...
78,...
74,...
84,...
109,...
115,...
116,...
120,...
12,...
53,...
51,...
59,...
117,...
121,...
118,...
124,...
20,...
61,...
56,...
67,...
119,...
125,...
116,...
120,...
22,...
90,...
88,...
96,...
107,...
111,...
118,...
124,...
27,...
99,...
93,...
105,...
109,...
115;...
29,...
38,...
31,...
39,...
10,...
33,...
28,...
37,...
35,...
47,...
36,...
48,...
19,...
42,...
34,...
46,...
50,...
58,...
51,...
59,...
12,...
53,...
49,...
57,...
55,...
66,...
56,...
67,...
20,...
61,...
54,...
65,...
69,...
76,...
70,...
77,...
6,...
71,...
68,...
75,...
73,...
83,...
74,...
84,...
14,...
78,...
72,...
82,...
50,...
58,...
31,...
39,...
9,...
32,...
68,...
75,...
55,...
66,...
36,...
48,...
18,...
41,...
72,...
82,...
86,...
95,...
88,...
96,...
22,...
90,...
85,...
94,...
92,...
104,...
93,...
105,...
27,...
99,...
91,...
103,...
107,...
111,...
70,...
77,...
12,...
53,...
106,...
110,...
109,...
115,...
74,...
84,...
20,...
61,...
108,...
114,...
117,...
121,...
51,...
59,...
8,...
52,...
116,...
120,...
119,...
125,...
56,...
67,...
13,...
60,...
118,...
124,...
107,...
111,...
88,...
96,...
21,...
89,...
116,...
120,...
109,...
115,...
93,...
105,...
26,...
98,...
118,...
124;...
30,...
11,...
37,...
34,...
38,...
35,...
39,...
36,...
40,...
1,...
46,...
43,...
47,...
44,...
48,...
45,...
32,...
18,...
57,...
54,...
58,...
55,...
59,...
56,...
41,...
15,...
65,...
62,...
66,...
63,...
67,...
64,...
33,...
19,...
75,...
72,...
76,...
73,...
77,...
74,...
42,...
16,...
82,...
79,...
83,...
80,...
84,...
81,...
53,...
20,...
75,...
72,...
58,...
55,...
39,...
36,...
61,...
17,...
82,...
79,...
66,...
63,...
48,...
45,...
87,...
23,...
94,...
91,...
95,...
92,...
96,...
93,...
97,...
3,...
103,...
100,...
104,...
101,...
105,...
102,...
89,...
26,...
110,...
108,...
111,...
109,...
77,...
74,...
98,...
24,...
114,...
112,...
115,...
113,...
84,...
81,...
90,...
27,...
120,...
118,...
121,...
119,...
59,...
56,...
99,...
25,...
124,...
122,...
125,...
123,...
67,...
64,...
53,...
20,...
120,...
118,...
111,...
109,...
96,...
93,...
61,...
17,...
124,...
122,...
115,...
113,...
105,...
102;...
37,...
34,...
32,...
18,...
39,...
36,...
38,...
35,...
46,...
43,...
41,...
15,...
48,...
45,...
47,...
44,...
57,...
54,...
52,...
13,...
59,...
56,...
58,...
55,...
65,...
62,...
60,...
4,...
67,...
64,...
66,...
63,...
75,...
72,...
53,...
20,...
77,...
74,...
76,...
73,...
82,...
79,...
61,...
17,...
84,...
81,...
83,...
80,...
75,...
72,...
33,...
19,...
39,...
36,...
58,...
55,...
82,...
79,...
42,...
16,...
48,...
45,...
66,...
63,...
94,...
91,...
89,...
26,...
96,...
93,...
95,...
92,...
103,...
100,...
98,...
24,...
105,...
102,...
104,...
101,...
110,...
108,...
71,...
14,...
77,...
74,...
111,...
109,...
114,...
112,...
78,...
2,...
84,...
81,...
115,...
113,...
120,...
118,...
53,...
20,...
59,...
56,...
121,...
119,...
124,...
122,...
61,...
17,...
67,...
64,...
125,...
123,...
120,...
118,...
90,...
27,...
96,...
93,...
111,...
109,...
124,...
122,...
99,...
25,...
105,...
102,...
115,...
113;...
38,...
35,...
39,...
36,...
33,...
19,...
37,...
34,...
47,...
44,...
48,...
45,...
42,...
16,...
46,...
43,...
58,...
55,...
59,...
56,...
53,...
20,...
57,...
54,...
66,...
63,...
67,...
64,...
61,...
17,...
65,...
62,...
76,...
73,...
77,...
74,...
71,...
14,...
75,...
72,...
83,...
80,...
84,...
81,...
78,...
2,...
82,...
79,...
58,...
55,...
39,...
36,...
32,...
18,...
75,...
72,...
66,...
63,...
48,...
45,...
41,...
15,...
82,...
79,...
95,...
92,...
96,...
93,...
90,...
27,...
94,...
91,...
104,...
101,...
105,...
102,...
99,...
25,...
103,...
100,...
111,...
109,...
77,...
74,...
53,...
20,...
110,...
108,...
115,...
113,...
84,...
81,...
61,...
17,...
114,...
112,...
121,...
119,...
59,...
56,...
52,...
13,...
120,...
118,...
125,...
123,...
67,...
64,...
60,...
4,...
124,...
122,...
111,...
109,...
96,...
93,...
89,...
26,...
120,...
118,...
115,...
113,...
105,...
102,...
98,...
24,...
124,...
122;...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
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
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
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
2,...
0,...
2,...
0,...
2,...
0,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
0,...
2,...
0,...
2,...
0,...
2,...
0,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
0,...
2,...
0,...
2,...
0,...
2,...
0,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
0,...
2,...
0,...
2,...
0,...
2,...
0,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
0,...
2,...
0,...
2,...
0,...
2,...
0,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
0,...
2,...
0,...
2,...
0,...
2,...
0,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
0,...
2,...
0,...
2,...
0,...
2,...
0,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
0,...
2,...
0,...
2,...
0,...
2,...
0,...
2,...
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
1,...
0,...
1,...
0,...
1,...
0,...
1,...
0,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
0,...
1,...
0,...
1,...
0,...
1,...
0,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
0,...
1,...
0,...
1,...
0,...
1,...
0,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
0,...
1,...
0,...
1,...
0,...
1,...
0,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
0,...
1,...
0,...
1,...
0,...
1,...
0,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
0,...
1,...
0,...
1,...
0,...
1,...
0,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
0,...
1,...
0,...
1,...
0,...
1,...
0,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
0,...
1,...
0,...
1,...
0,...
1,...
0;...
0,...
0,...
0,...
0,...
3,...
3,...
3,...
3,...
0,...
0,...
0,...
0,...
3,...
3,...
3,...
3,...
0,...
0,...
0,...
0,...
3,...
3,...
3,...
3,...
0,...
0,...
0,...
0,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
0,...
0,...
0,...
0,...
3,...
3,...
3,...
3,...
0,...
0,...
0,...
0,...
3,...
3,...
3,...
3,...
0,...
0,...
0,...
0,...
3,...
3,...
3,...
3,...
0,...
0,...
0,...
0,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3,...
3;...
0,...
0,...
4,...
4,...
0,...
0,...
4,...
4,...
0,...
0,...
4,...
4,...
0,...
0,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
0,...
0,...
4,...
4,...
0,...
0,...
4,...
4,...
0,...
0,...
4,...
4,...
0,...
0,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
0,...
0,...
4,...
4,...
0,...
0,...
4,...
4,...
0,...
0,...
4,...
4,...
0,...
0,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
0,...
0,...
4,...
4,...
0,...
0,...
4,...
4,...
0,...
0,...
4,...
4,...
0,...
0,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4;...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5,...
5;...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
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
1,...
0,...
3,...
0,...
5,...
0,...
7,...
2,...
9,...
4,...
11,...
6,...
13,...
8,...
15,...
0,...
17,...
0,...
19,...
0,...
21,...
0,...
23,...
18,...
25,...
20,...
27,...
22,...
29,...
24,...
31,...
0,...
33,...
0,...
35,...
0,...
37,...
0,...
39,...
34,...
41,...
36,...
43,...
38,...
45,...
40,...
47,...
0,...
49,...
0,...
51,...
0,...
53,...
0,...
55,...
50,...
57,...
52,...
59,...
54,...
61,...
56,...
63,...
0,...
65,...
0,...
67,...
0,...
69,...
0,...
71,...
66,...
73,...
68,...
75,...
70,...
77,...
72,...
79,...
0,...
81,...
0,...
83,...
0,...
85,...
0,...
87,...
82,...
89,...
84,...
91,...
86,...
93,...
88,...
95,...
0,...
97,...
0,...
99,...
0,...
101,...
0,...
103,...
98,...
105,...
100,...
107,...
102,...
109,...
104,...
111,...
0,...
113,...
0,...
115,...
0,...
117,...
0,...
119,...
114,...
121,...
116,...
123,...
118,...
125,...
120,...
127;...
2,...
9,...
4,...
11,...
6,...
13,...
8,...
15,...
10,...
0,...
12,...
0,...
14,...
0,...
16,...
0,...
18,...
25,...
20,...
27,...
22,...
29,...
24,...
31,...
26,...
0,...
28,...
0,...
30,...
0,...
32,...
0,...
34,...
41,...
36,...
43,...
38,...
45,...
40,...
47,...
42,...
0,...
44,...
0,...
46,...
0,...
48,...
0,...
50,...
57,...
52,...
59,...
54,...
61,...
56,...
63,...
58,...
0,...
60,...
0,...
62,...
0,...
64,...
0,...
66,...
73,...
68,...
75,...
70,...
77,...
72,...
79,...
74,...
0,...
76,...
0,...
78,...
0,...
80,...
0,...
82,...
89,...
84,...
91,...
86,...
93,...
88,...
95,...
90,...
0,...
92,...
0,...
94,...
0,...
96,...
0,...
98,...
105,...
100,...
107,...
102,...
109,...
104,...
111,...
106,...
0,...
108,...
0,...
110,...
0,...
112,...
0,...
114,...
121,...
116,...
123,...
118,...
125,...
120,...
127,...
122,...
0,...
124,...
0,...
126,...
0,...
128,...
0;...
0,...
0,...
0,...
0,...
7,...
8,...
5,...
6,...
0,...
0,...
0,...
0,...
15,...
16,...
13,...
14,...
0,...
0,...
0,...
0,...
23,...
24,...
21,...
22,...
0,...
0,...
0,...
0,...
31,...
32,...
29,...
30,...
51,...
52,...
49,...
50,...
39,...
40,...
37,...
38,...
59,...
60,...
57,...
58,...
47,...
48,...
45,...
46,...
35,...
36,...
33,...
34,...
55,...
56,...
53,...
54,...
43,...
44,...
41,...
42,...
63,...
64,...
61,...
62,...
0,...
0,...
0,...
0,...
71,...
72,...
69,...
70,...
0,...
0,...
0,...
0,...
79,...
80,...
77,...
78,...
0,...
0,...
0,...
0,...
87,...
88,...
85,...
86,...
0,...
0,...
0,...
0,...
95,...
96,...
93,...
94,...
115,...
116,...
113,...
114,...
103,...
104,...
101,...
102,...
123,...
124,...
121,...
122,...
111,...
112,...
109,...
110,...
99,...
100,...
97,...
98,...
119,...
120,...
117,...
118,...
107,...
108,...
105,...
106,...
127,...
128,...
125,...
126;...
0,...
0,...
7,...
8,...
0,...
0,...
3,...
4,...
0,...
0,...
15,...
16,...
0,...
0,...
11,...
12,...
53,...
54,...
23,...
24,...
49,...
50,...
19,...
20,...
61,...
62,...
31,...
32,...
57,...
58,...
27,...
28,...
0,...
0,...
39,...
40,...
0,...
0,...
35,...
36,...
0,...
0,...
47,...
48,...
0,...
0,...
43,...
44,...
21,...
22,...
55,...
56,...
17,...
18,...
51,...
52,...
29,...
30,...
63,...
64,...
25,...
26,...
59,...
60,...
0,...
0,...
71,...
72,...
0,...
0,...
67,...
68,...
0,...
0,...
79,...
80,...
0,...
0,...
75,...
76,...
117,...
118,...
87,...
88,...
113,...
114,...
83,...
84,...
125,...
126,...
95,...
96,...
121,...
122,...
91,...
92,...
0,...
0,...
103,...
104,...
0,...
0,...
99,...
100,...
0,...
0,...
111,...
112,...
0,...
0,...
107,...
108,...
85,...
86,...
119,...
120,...
81,...
82,...
115,...
116,...
93,...
94,...
127,...
128,...
89,...
90,...
123,...
124;...
7,...
8,...
53,...
54,...
51,...
52,...
1,...
2,...
15,...
16,...
61,...
62,...
59,...
60,...
9,...
10,...
23,...
24,...
101,...
102,...
99,...
100,...
17,...
18,...
31,...
32,...
109,...
110,...
107,...
108,...
25,...
26,...
39,...
40,...
85,...
86,...
83,...
84,...
33,...
34,...
47,...
48,...
93,...
94,...
91,...
92,...
41,...
42,...
55,...
56,...
5,...
6,...
3,...
4,...
49,...
50,...
63,...
64,...
13,...
14,...
11,...
12,...
57,...
58,...
71,...
72,...
117,...
118,...
115,...
116,...
65,...
66,...
79,...
80,...
125,...
126,...
123,...
124,...
73,...
74,...
87,...
88,...
37,...
38,...
35,...
36,...
81,...
82,...
95,...
96,...
45,...
46,...
43,...
44,...
89,...
90,...
103,...
104,...
21,...
22,...
19,...
20,...
97,...
98,...
111,...
112,...
29,...
30,...
27,...
28,...
105,...
106,...
119,...
120,...
69,...
70,...
67,...
68,...
113,...
114,...
127,...
128,...
77,...
78,...
75,...
76,...
121,...
122;...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
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
