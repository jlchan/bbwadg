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
-0.5; ...
0; ...
0.5; ...
1; ...
1; ...
1; ...
1; ...
-1; ...
0.5; ...
0; ...
-0.5; ...
-1; ...
-1; ...
-1; ...
-0.588646; ...
0.211667; ...
0.588648; ...
0.626762; ...
-0.626772; ...
0.188415; ...
-0.626379; ...
0.626387; ...
-0.188414; ...
0.343943; ...
-0.211677; ...
0.718062; ...
-0.718064; ...
0.15189; ...
-0.151889; ...
-0.343943; ...
-1; ...
1; ...
-0.5; ...
0; ...
0.5; ...
1; ...
1; ...
1; ...
1; ...
-1; ...
0.5; ...
0; ...
-0.5; ...
-1; ...
-1; ...
-1; ...
-0.588646; ...
0.211667; ...
0.588648; ...
0.626762; ...
-0.626772; ...
0.188415; ...
-0.626379; ...
0.626387; ...
-0.188414; ...
0.343943; ...
-0.211677; ...
0.718062; ...
-0.718064; ...
0.15189; ...
-0.151889; ...
-0.343943; ...
-1; ...
1; ...
-0.5; ...
0; ...
0.5; ...
1; ...
1; ...
1; ...
1; ...
-1; ...
0.5; ...
0; ...
-0.5; ...
-1; ...
-1; ...
-1; ...
-0.588646; ...
0.211667; ...
0.588648; ...
0.626762; ...
-0.626772; ...
0.188415; ...
-0.626379; ...
0.626387; ...
-0.188414; ...
0.343943; ...
-0.211677; ...
0.718062; ...
-0.718064; ...
0.15189; ...
-0.151889; ...
-0.343943; ...
-1; ...
1; ...
-0.5; ...
0; ...
0.5; ...
1; ...
1; ...
1; ...
1; ...
-1; ...
0.5; ...
0; ...
-0.5; ...
-1; ...
-1; ...
-1; ...
-0.588646; ...
0.211667; ...
0.588648; ...
0.626762; ...
-0.626772; ...
0.188415; ...
-0.626379; ...
0.626387; ...
-0.188414; ...
0.343943; ...
-0.211677; ...
0.718062; ...
-0.718064; ...
0.15189; ...
-0.151889; ...
-0.343943; ...
-1; ...
1; ...
-0.5; ...
0; ...
0.5; ...
1; ...
1; ...
1; ...
1; ...
-1; ...
0.5; ...
0; ...
-0.5; ...
-1; ...
-1; ...
-1; ...
-0.588646; ...
0.211667; ...
0.588648; ...
0.626762; ...
-0.626772; ...
0.188415; ...
-0.626379; ...
0.626387; ...
-0.188414; ...
0.343943; ...
-0.211677; ...
0.718062; ...
-0.718064; ...
0.15189; ...
-0.151889; ...
-0.343943; ...
];
VY = [ ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
1; ...
-0.5; ...
0; ...
0.5; ...
1; ...
1; ...
1; ...
1; ...
0.5; ...
0; ...
-0.5; ...
-0.598992; ...
0.685733; ...
0.599003; ...
-0.606914; ...
0.606902; ...
-0.620334; ...
-0.163523; ...
0.163515; ...
0.62033; ...
-0.209402; ...
-0.685741; ...
-0.230386; ...
0.230388; ...
0.248404; ...
-0.2484; ...
0.209398; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
1; ...
-0.5; ...
0; ...
0.5; ...
1; ...
1; ...
1; ...
1; ...
0.5; ...
0; ...
-0.5; ...
-0.598992; ...
0.685733; ...
0.599003; ...
-0.606914; ...
0.606902; ...
-0.620334; ...
-0.163523; ...
0.163515; ...
0.62033; ...
-0.209402; ...
-0.685741; ...
-0.230386; ...
0.230388; ...
0.248404; ...
-0.2484; ...
0.209398; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
1; ...
-0.5; ...
0; ...
0.5; ...
1; ...
1; ...
1; ...
1; ...
0.5; ...
0; ...
-0.5; ...
-0.598992; ...
0.685733; ...
0.599003; ...
-0.606914; ...
0.606902; ...
-0.620334; ...
-0.163523; ...
0.163515; ...
0.62033; ...
-0.209402; ...
-0.685741; ...
-0.230386; ...
0.230388; ...
0.248404; ...
-0.2484; ...
0.209398; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
1; ...
-0.5; ...
0; ...
0.5; ...
1; ...
1; ...
1; ...
1; ...
0.5; ...
0; ...
-0.5; ...
-0.598992; ...
0.685733; ...
0.599003; ...
-0.606914; ...
0.606902; ...
-0.620334; ...
-0.163523; ...
0.163515; ...
0.62033; ...
-0.209402; ...
-0.685741; ...
-0.230386; ...
0.230388; ...
0.248404; ...
-0.2484; ...
0.209398; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
1; ...
-0.5; ...
0; ...
0.5; ...
1; ...
1; ...
1; ...
1; ...
0.5; ...
0; ...
-0.5; ...
-0.598992; ...
0.685733; ...
0.599003; ...
-0.606914; ...
0.606902; ...
-0.620334; ...
-0.163523; ...
0.163515; ...
0.62033; ...
-0.209402; ...
-0.685741; ...
-0.230386; ...
0.230388; ...
0.248404; ...
-0.2484; ...
0.209398; ...
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
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-1; ...
-0.625; ...
-0.375; ...
-0.625; ...
-0.375; ...
-0.625; ...
-0.375; ...
-0.625; ...
-0.375; ...
-0.625; ...
-0.375; ...
-0.625; ...
-0.375; ...
-0.625; ...
-0.375; ...
-0.625; ...
-0.375; ...
-0.625; ...
-0.375; ...
-0.625; ...
-0.375; ...
-0.625; ...
-0.375; ...
-0.625; ...
-0.375; ...
-0.625; ...
-0.375; ...
-0.625; ...
-0.375; ...
-0.625; ...
-0.375; ...
-0.625; ...
-0.375; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
-0; ...
0.375; ...
0.625; ...
0.375; ...
0.625; ...
0.375; ...
0.625; ...
0.375; ...
0.625; ...
0.375; ...
0.625; ...
0.375; ...
0.625; ...
0.375; ...
0.625; ...
0.375; ...
0.625; ...
0.375; ...
0.625; ...
0.375; ...
0.625; ...
0.375; ...
0.625; ...
0.375; ...
0.625; ...
0.375; ...
0.625; ...
0.375; ...
0.625; ...
0.375; ...
0.625; ...
0.375; ...
0.625; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
1; ...
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
16,...
17,...
9,...
19,...
5,...
20,...
13,...
21,...
5,...
5,...
16,...
16,...
9,...
9,...
13,...
13,...
19,...
18,...
17,...
27,...
20,...
28,...
21,...
29,...
18,...
27,...
28,...
29,...
20,...
20,...
24,...
30,...
23,...
31,...
21,...
21,...
28,...
24,...
29,...
23,...
27,...
22,...
18,...
25,...
30,...
31,...
48,...
49,...
41,...
51,...
37,...
52,...
45,...
53,...
37,...
37,...
48,...
48,...
41,...
41,...
45,...
45,...
51,...
50,...
49,...
59,...
52,...
60,...
53,...
61,...
50,...
59,...
60,...
61,...
52,...
52,...
56,...
62,...
55,...
63,...
53,...
53,...
60,...
56,...
61,...
55,...
59,...
54,...
50,...
57,...
62,...
63,...
80,...
81,...
73,...
83,...
69,...
84,...
77,...
85,...
69,...
69,...
80,...
80,...
73,...
73,...
77,...
77,...
83,...
82,...
81,...
91,...
84,...
92,...
85,...
93,...
82,...
91,...
92,...
93,...
84,...
84,...
88,...
94,...
87,...
95,...
85,...
85,...
92,...
88,...
93,...
87,...
91,...
86,...
82,...
89,...
94,...
95,...
112,...
113,...
105,...
115,...
101,...
116,...
109,...
117,...
101,...
101,...
112,...
112,...
105,...
105,...
109,...
109,...
115,...
114,...
113,...
123,...
116,...
124,...
117,...
125,...
114,...
123,...
124,...
125,...
116,...
116,...
120,...
126,...
119,...
127,...
117,...
117,...
124,...
120,...
125,...
119,...
123,...
118,...
114,...
121,...
126,...
127;...
1,...
1,...
6,...
6,...
2,...
2,...
10,...
10,...
20,...
22,...
17,...
23,...
19,...
24,...
21,...
25,...
11,...
11,...
3,...
3,...
7,...
7,...
14,...
14,...
12,...
4,...
8,...
15,...
28,...
26,...
19,...
19,...
17,...
17,...
29,...
32,...
24,...
30,...
23,...
31,...
22,...
26,...
25,...
32,...
32,...
26,...
33,...
33,...
38,...
38,...
34,...
34,...
42,...
42,...
52,...
54,...
49,...
55,...
51,...
56,...
53,...
57,...
43,...
43,...
35,...
35,...
39,...
39,...
46,...
46,...
44,...
36,...
40,...
47,...
60,...
58,...
51,...
51,...
49,...
49,...
61,...
64,...
56,...
62,...
55,...
63,...
54,...
58,...
57,...
64,...
64,...
58,...
65,...
65,...
70,...
70,...
66,...
66,...
74,...
74,...
84,...
86,...
81,...
87,...
83,...
88,...
85,...
89,...
75,...
75,...
67,...
67,...
71,...
71,...
78,...
78,...
76,...
68,...
72,...
79,...
92,...
90,...
83,...
83,...
81,...
81,...
93,...
96,...
88,...
94,...
87,...
95,...
86,...
90,...
89,...
96,...
96,...
90,...
97,...
97,...
102,...
102,...
98,...
98,...
106,...
106,...
116,...
118,...
113,...
119,...
115,...
120,...
117,...
121,...
107,...
107,...
99,...
99,...
103,...
103,...
110,...
110,...
108,...
100,...
104,...
111,...
124,...
122,...
115,...
115,...
113,...
113,...
125,...
128,...
120,...
126,...
119,...
127,...
118,...
122,...
121,...
128,...
128,...
122;...
17,...
3,...
19,...
11,...
20,...
7,...
21,...
14,...
22,...
4,...
23,...
15,...
24,...
8,...
25,...
12,...
18,...
12,...
27,...
4,...
28,...
8,...
29,...
15,...
25,...
22,...
24,...
23,...
26,...
22,...
30,...
18,...
31,...
27,...
32,...
25,...
26,...
26,...
32,...
32,...
31,...
31,...
30,...
30,...
31,...
30,...
49,...
35,...
51,...
43,...
52,...
39,...
53,...
46,...
54,...
36,...
55,...
47,...
56,...
40,...
57,...
44,...
50,...
44,...
59,...
36,...
60,...
40,...
61,...
47,...
57,...
54,...
56,...
55,...
58,...
54,...
62,...
50,...
63,...
59,...
64,...
57,...
58,...
58,...
64,...
64,...
63,...
63,...
62,...
62,...
63,...
62,...
81,...
67,...
83,...
75,...
84,...
71,...
85,...
78,...
86,...
68,...
87,...
79,...
88,...
72,...
89,...
76,...
82,...
76,...
91,...
68,...
92,...
72,...
93,...
79,...
89,...
86,...
88,...
87,...
90,...
86,...
94,...
82,...
95,...
91,...
96,...
89,...
90,...
90,...
96,...
96,...
95,...
95,...
94,...
94,...
95,...
94,...
113,...
99,...
115,...
107,...
116,...
103,...
117,...
110,...
118,...
100,...
119,...
111,...
120,...
104,...
121,...
108,...
114,...
108,...
123,...
100,...
124,...
104,...
125,...
111,...
121,...
118,...
120,...
119,...
122,...
118,...
126,...
114,...
127,...
123,...
128,...
121,...
122,...
122,...
128,...
128,...
127,...
127,...
126,...
126,...
127,...
126;...
48,...
49,...
41,...
51,...
37,...
52,...
45,...
53,...
37,...
37,...
48,...
48,...
41,...
41,...
45,...
45,...
51,...
50,...
49,...
59,...
52,...
60,...
53,...
61,...
50,...
59,...
60,...
61,...
52,...
52,...
56,...
62,...
55,...
63,...
53,...
53,...
60,...
56,...
61,...
55,...
59,...
54,...
50,...
57,...
62,...
63,...
80,...
81,...
73,...
83,...
69,...
84,...
77,...
85,...
69,...
69,...
80,...
80,...
73,...
73,...
77,...
77,...
83,...
82,...
81,...
91,...
84,...
92,...
85,...
93,...
82,...
91,...
92,...
93,...
84,...
84,...
88,...
94,...
87,...
95,...
85,...
85,...
92,...
88,...
93,...
87,...
91,...
86,...
82,...
89,...
94,...
95,...
112,...
113,...
105,...
115,...
101,...
116,...
109,...
117,...
101,...
101,...
112,...
112,...
105,...
105,...
109,...
109,...
115,...
114,...
113,...
123,...
116,...
124,...
117,...
125,...
114,...
123,...
124,...
125,...
116,...
116,...
120,...
126,...
119,...
127,...
117,...
117,...
124,...
120,...
125,...
119,...
123,...
118,...
114,...
121,...
126,...
127,...
144,...
145,...
137,...
147,...
133,...
148,...
141,...
149,...
133,...
133,...
144,...
144,...
137,...
137,...
141,...
141,...
147,...
146,...
145,...
155,...
148,...
156,...
149,...
157,...
146,...
155,...
156,...
157,...
148,...
148,...
152,...
158,...
151,...
159,...
149,...
149,...
156,...
152,...
157,...
151,...
155,...
150,...
146,...
153,...
158,...
159;...
33,...
33,...
38,...
38,...
34,...
34,...
42,...
42,...
52,...
54,...
49,...
55,...
51,...
56,...
53,...
57,...
43,...
43,...
35,...
35,...
39,...
39,...
46,...
46,...
44,...
36,...
40,...
47,...
60,...
58,...
51,...
51,...
49,...
49,...
61,...
64,...
56,...
62,...
55,...
63,...
54,...
58,...
57,...
64,...
64,...
58,...
65,...
65,...
70,...
70,...
66,...
66,...
74,...
74,...
84,...
86,...
81,...
87,...
83,...
88,...
85,...
89,...
75,...
75,...
67,...
67,...
71,...
71,...
78,...
78,...
76,...
68,...
72,...
79,...
92,...
90,...
83,...
83,...
81,...
81,...
93,...
96,...
88,...
94,...
87,...
95,...
86,...
90,...
89,...
96,...
96,...
90,...
97,...
97,...
102,...
102,...
98,...
98,...
106,...
106,...
116,...
118,...
113,...
119,...
115,...
120,...
117,...
121,...
107,...
107,...
99,...
99,...
103,...
103,...
110,...
110,...
108,...
100,...
104,...
111,...
124,...
122,...
115,...
115,...
113,...
113,...
125,...
128,...
120,...
126,...
119,...
127,...
118,...
122,...
121,...
128,...
128,...
122,...
129,...
129,...
134,...
134,...
130,...
130,...
138,...
138,...
148,...
150,...
145,...
151,...
147,...
152,...
149,...
153,...
139,...
139,...
131,...
131,...
135,...
135,...
142,...
142,...
140,...
132,...
136,...
143,...
156,...
154,...
147,...
147,...
145,...
145,...
157,...
160,...
152,...
158,...
151,...
159,...
150,...
154,...
153,...
160,...
160,...
154;...
49,...
35,...
51,...
43,...
52,...
39,...
53,...
46,...
54,...
36,...
55,...
47,...
56,...
40,...
57,...
44,...
50,...
44,...
59,...
36,...
60,...
40,...
61,...
47,...
57,...
54,...
56,...
55,...
58,...
54,...
62,...
50,...
63,...
59,...
64,...
57,...
58,...
58,...
64,...
64,...
63,...
63,...
62,...
62,...
63,...
62,...
81,...
67,...
83,...
75,...
84,...
71,...
85,...
78,...
86,...
68,...
87,...
79,...
88,...
72,...
89,...
76,...
82,...
76,...
91,...
68,...
92,...
72,...
93,...
79,...
89,...
86,...
88,...
87,...
90,...
86,...
94,...
82,...
95,...
91,...
96,...
89,...
90,...
90,...
96,...
96,...
95,...
95,...
94,...
94,...
95,...
94,...
113,...
99,...
115,...
107,...
116,...
103,...
117,...
110,...
118,...
100,...
119,...
111,...
120,...
104,...
121,...
108,...
114,...
108,...
123,...
100,...
124,...
104,...
125,...
111,...
121,...
118,...
120,...
119,...
122,...
118,...
126,...
114,...
127,...
123,...
128,...
121,...
122,...
122,...
128,...
128,...
127,...
127,...
126,...
126,...
127,...
126,...
145,...
131,...
147,...
139,...
148,...
135,...
149,...
142,...
150,...
132,...
151,...
143,...
152,...
136,...
153,...
140,...
146,...
140,...
155,...
132,...
156,...
136,...
157,...
143,...
153,...
150,...
152,...
151,...
154,...
150,...
158,...
146,...
159,...
155,...
160,...
153,...
154,...
154,...
160,...
160,...
159,...
159,...
158,...
158,...
159,...
158;...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
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
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
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
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
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
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
2,...
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
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
1,...
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
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
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
5,...
0,...
5,...
0,...
5,...
0,...
5,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
5,...
4,...
5,...
4,...
5,...
4,...
5,...
4,...
4,...
4,...
4,...
4,...
4,...
5,...
5,...
5,...
5,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
5,...
4,...
5,...
5,...
5,...
0,...
5,...
0,...
5,...
0,...
5,...
0,...
5,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
5,...
4,...
5,...
4,...
5,...
4,...
5,...
4,...
4,...
4,...
4,...
4,...
4,...
5,...
5,...
5,...
5,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
5,...
4,...
5,...
5,...
5,...
0,...
5,...
0,...
5,...
0,...
5,...
0,...
5,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
5,...
4,...
5,...
4,...
5,...
4,...
5,...
4,...
4,...
4,...
4,...
4,...
4,...
5,...
5,...
5,...
5,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
5,...
4,...
5,...
5,...
5,...
0,...
5,...
0,...
5,...
0,...
5,...
0,...
5,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
5,...
4,...
5,...
4,...
5,...
4,...
5,...
4,...
4,...
4,...
4,...
4,...
4,...
5,...
5,...
5,...
5,...
4,...
4,...
4,...
4,...
4,...
4,...
4,...
5,...
4,...
5,...
5,...
5;...
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
3,...
0,...
3,...
0,...
3,...
0,...
5,...
3,...
5,...
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
5,...
3,...
4,...
3,...
4,...
3,...
5,...
5,...
5,...
5,...
5,...
4,...
5,...
4,...
5,...
4,...
4,...
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
3,...
0,...
3,...
0,...
3,...
0,...
5,...
3,...
5,...
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
5,...
3,...
4,...
3,...
4,...
3,...
5,...
5,...
5,...
5,...
5,...
4,...
5,...
4,...
5,...
4,...
4,...
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
3,...
0,...
3,...
0,...
3,...
0,...
5,...
3,...
5,...
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
5,...
3,...
4,...
3,...
4,...
3,...
5,...
5,...
5,...
5,...
5,...
4,...
5,...
4,...
5,...
4,...
4,...
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
3,...
0,...
3,...
0,...
3,...
0,...
5,...
3,...
5,...
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
5,...
3,...
4,...
3,...
4,...
3,...
5,...
5,...
5,...
5,...
5,...
4,...
5,...
4,...
5,...
4,...
4;...
3,...
0,...
3,...
0,...
3,...
0,...
3,...
0,...
4,...
5,...
3,...
5,...
3,...
5,...
4,...
5,...
3,...
0,...
3,...
0,...
3,...
0,...
3,...
0,...
5,...
5,...
5,...
5,...
4,...
3,...
3,...
4,...
3,...
4,...
4,...
3,...
4,...
5,...
4,...
5,...
4,...
3,...
4,...
3,...
5,...
5,...
3,...
0,...
3,...
0,...
3,...
0,...
3,...
0,...
4,...
5,...
3,...
5,...
3,...
5,...
4,...
5,...
3,...
0,...
3,...
0,...
3,...
0,...
3,...
0,...
5,...
5,...
5,...
5,...
4,...
3,...
3,...
4,...
3,...
4,...
4,...
3,...
4,...
5,...
4,...
5,...
4,...
3,...
4,...
3,...
5,...
5,...
3,...
0,...
3,...
0,...
3,...
0,...
3,...
0,...
4,...
5,...
3,...
5,...
3,...
5,...
4,...
5,...
3,...
0,...
3,...
0,...
3,...
0,...
3,...
0,...
5,...
5,...
5,...
5,...
4,...
3,...
3,...
4,...
3,...
4,...
4,...
3,...
4,...
5,...
4,...
5,...
4,...
3,...
4,...
3,...
5,...
5,...
3,...
0,...
3,...
0,...
3,...
0,...
3,...
0,...
4,...
5,...
3,...
5,...
3,...
5,...
4,...
5,...
3,...
0,...
3,...
0,...
3,...
0,...
3,...
0,...
5,...
5,...
5,...
5,...
4,...
3,...
3,...
4,...
3,...
4,...
4,...
3,...
4,...
5,...
4,...
5,...
4,...
3,...
4,...
3,...
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
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
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
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
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
8,...
9,...
10,...
11,...
12,...
13,...
14,...
15,...
16,...
17,...
18,...
19,...
20,...
21,...
22,...
23,...
24,...
25,...
26,...
27,...
28,...
29,...
30,...
31,...
32,...
33,...
34,...
35,...
36,...
37,...
38,...
39,...
40,...
41,...
42,...
43,...
44,...
45,...
46,...
47,...
48,...
49,...
50,...
51,...
52,...
53,...
54,...
55,...
56,...
57,...
58,...
59,...
60,...
61,...
62,...
63,...
64,...
65,...
66,...
67,...
68,...
69,...
70,...
71,...
72,...
73,...
74,...
75,...
76,...
77,...
78,...
79,...
80,...
81,...
82,...
83,...
84,...
85,...
86,...
87,...
88,...
89,...
90,...
91,...
92,...
93,...
94,...
95,...
96,...
97,...
98,...
99,...
100,...
101,...
102,...
103,...
104,...
105,...
106,...
107,...
108,...
109,...
110,...
111,...
112,...
113,...
114,...
115,...
116,...
117,...
118,...
119,...
120,...
121,...
122,...
123,...
124,...
125,...
126,...
127,...
128,...
129,...
130,...
131,...
132,...
133,...
134,...
135,...
136,...
137,...
138;...
47,...
48,...
49,...
50,...
51,...
52,...
53,...
54,...
55,...
56,...
57,...
58,...
59,...
60,...
61,...
62,...
63,...
64,...
65,...
66,...
67,...
68,...
69,...
70,...
71,...
72,...
73,...
74,...
75,...
76,...
77,...
78,...
79,...
80,...
81,...
82,...
83,...
84,...
85,...
86,...
87,...
88,...
89,...
90,...
91,...
92,...
93,...
94,...
95,...
96,...
97,...
98,...
99,...
100,...
101,...
102,...
103,...
104,...
105,...
106,...
107,...
108,...
109,...
110,...
111,...
112,...
113,...
114,...
115,...
116,...
117,...
118,...
119,...
120,...
121,...
122,...
123,...
124,...
125,...
126,...
127,...
128,...
129,...
130,...
131,...
132,...
133,...
134,...
135,...
136,...
137,...
138,...
139,...
140,...
141,...
142,...
143,...
144,...
145,...
146,...
147,...
148,...
149,...
150,...
151,...
152,...
153,...
154,...
155,...
156,...
157,...
158,...
159,...
160,...
161,...
162,...
163,...
164,...
165,...
166,...
167,...
168,...
169,...
170,...
171,...
172,...
173,...
174,...
175,...
176,...
177,...
178,...
179,...
180,...
181,...
182,...
183,...
184,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
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
1,...
0,...
3,...
0,...
5,...
0,...
7,...
5,...
9,...
1,...
11,...
3,...
13,...
7,...
15,...
4,...
17,...
2,...
19,...
6,...
21,...
8,...
23,...
18,...
20,...
22,...
24,...
21,...
29,...
13,...
31,...
11,...
33,...
23,...
35,...
27,...
31,...
28,...
33,...
26,...
30,...
25,...
36,...
44,...
42,...
0,...
47,...
0,...
49,...
0,...
51,...
0,...
53,...
51,...
55,...
47,...
57,...
49,...
59,...
53,...
61,...
50,...
63,...
48,...
65,...
52,...
67,...
54,...
69,...
64,...
66,...
68,...
70,...
67,...
75,...
59,...
77,...
57,...
79,...
69,...
81,...
73,...
77,...
74,...
79,...
72,...
76,...
71,...
82,...
90,...
88,...
0,...
93,...
0,...
95,...
0,...
97,...
0,...
99,...
97,...
101,...
93,...
103,...
95,...
105,...
99,...
107,...
96,...
109,...
94,...
111,...
98,...
113,...
100,...
115,...
110,...
112,...
114,...
116,...
113,...
121,...
105,...
123,...
103,...
125,...
115,...
127,...
119,...
123,...
120,...
125,...
118,...
122,...
117,...
128,...
136,...
134,...
0,...
139,...
0,...
141,...
0,...
143,...
0,...
145,...
143,...
147,...
139,...
149,...
141,...
151,...
145,...
153,...
142,...
155,...
140,...
157,...
144,...
159,...
146,...
161,...
156,...
158,...
160,...
162,...
159,...
167,...
151,...
169,...
149,...
171,...
161,...
173,...
165,...
169,...
166,...
171,...
164,...
168,...
163,...
174,...
182,...
180;...
11,...
19,...
13,...
17,...
9,...
21,...
15,...
23,...
10,...
0,...
12,...
0,...
14,...
0,...
16,...
0,...
32,...
25,...
34,...
26,...
29,...
27,...
35,...
28,...
43,...
41,...
37,...
39,...
30,...
9,...
38,...
43,...
40,...
41,...
36,...
15,...
29,...
37,...
35,...
39,...
34,...
41,...
32,...
43,...
46,...
45,...
57,...
65,...
59,...
63,...
55,...
67,...
61,...
69,...
56,...
0,...
58,...
0,...
60,...
0,...
62,...
0,...
78,...
71,...
80,...
72,...
75,...
73,...
81,...
74,...
89,...
87,...
83,...
85,...
76,...
55,...
84,...
89,...
86,...
87,...
82,...
61,...
75,...
83,...
81,...
85,...
80,...
87,...
78,...
89,...
92,...
91,...
103,...
111,...
105,...
109,...
101,...
113,...
107,...
115,...
102,...
0,...
104,...
0,...
106,...
0,...
108,...
0,...
124,...
117,...
126,...
118,...
121,...
119,...
127,...
120,...
135,...
133,...
129,...
131,...
122,...
101,...
130,...
135,...
132,...
133,...
128,...
107,...
121,...
129,...
127,...
131,...
126,...
133,...
124,...
135,...
138,...
137,...
149,...
157,...
151,...
155,...
147,...
159,...
153,...
161,...
148,...
0,...
150,...
0,...
152,...
0,...
154,...
0,...
170,...
163,...
172,...
164,...
167,...
165,...
173,...
166,...
181,...
179,...
175,...
177,...
168,...
147,...
176,...
181,...
178,...
179,...
174,...
153,...
167,...
175,...
173,...
177,...
172,...
179,...
170,...
181,...
184,...
183;...
2,...
0,...
4,...
0,...
6,...
0,...
8,...
0,...
30,...
26,...
33,...
28,...
31,...
27,...
36,...
25,...
18,...
0,...
20,...
0,...
22,...
0,...
24,...
0,...
16,...
10,...
14,...
12,...
37,...
42,...
32,...
17,...
34,...
19,...
39,...
44,...
38,...
46,...
40,...
45,...
42,...
46,...
44,...
45,...
40,...
38,...
48,...
0,...
50,...
0,...
52,...
0,...
54,...
0,...
76,...
72,...
79,...
74,...
77,...
73,...
82,...
71,...
64,...
0,...
66,...
0,...
68,...
0,...
70,...
0,...
62,...
56,...
60,...
58,...
83,...
88,...
78,...
63,...
80,...
65,...
85,...
90,...
84,...
92,...
86,...
91,...
88,...
92,...
90,...
91,...
86,...
84,...
94,...
0,...
96,...
0,...
98,...
0,...
100,...
0,...
122,...
118,...
125,...
120,...
123,...
119,...
128,...
117,...
110,...
0,...
112,...
0,...
114,...
0,...
116,...
0,...
108,...
102,...
106,...
104,...
129,...
134,...
124,...
109,...
126,...
111,...
131,...
136,...
130,...
138,...
132,...
137,...
134,...
138,...
136,...
137,...
132,...
130,...
140,...
0,...
142,...
0,...
144,...
0,...
146,...
0,...
168,...
164,...
171,...
166,...
169,...
165,...
174,...
163,...
156,...
0,...
158,...
0,...
160,...
0,...
162,...
0,...
154,...
148,...
152,...
150,...
175,...
180,...
170,...
155,...
172,...
157,...
177,...
182,...
176,...
184,...
178,...
183,...
180,...
184,...
182,...
183,...
178,...
176;...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
0,...
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
