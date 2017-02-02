// Inclusion radius (i.e. side length is 2r)
rx = .250000;
ry = .250000;
rz = .250000;


// Inclusion is centered at (0,0,z0)
x0 = 0;
y0 = 0;
z0 = -.00000;

// PML thickness
d = 0.00000;

// Element size
h = 0.125;

//====================================
// Cube 1
//====================================
Point(1) = {-rx+x0, -ry+y0, -rz+z0, h};
Point(2) = {rx+x0, -ry+y0, -rz+z0, h};
Point(3) = {-rx+x0, ry+y0, -rz+z0, h};
Point(4) = {rx+x0, ry+y0, -rz+z0, h};
Point(5) = {rx+x0, ry+y0, rz+z0, h};
Point(6) = {rx+x0, -ry+y0, rz+z0, h};
Point(7) = {-rx+x0, ry+y0, rz+z0, h};
Point(8) = {-rx+x0, -ry+y0, rz+z0, h};

Line(1) = {3, 7};
Line(2) = {7, 5};
Line(3) = {5, 4};
Line(4) = {4, 3};
Line(5) = {3, 1};
Line(6) = {2, 4};
Line(7) = {2, 6};
Line(8) = {6, 8};
Line(9) = {8, 1};
Line(10) = {1, 2};
Line(11) = {8, 7};
Line(12) = {6, 5};

Line Loop(13) = {7, 8, 9, 10};
Plane Surface(14) = {13};
Line Loop(15) = {6, 4, 5, 10};
Plane Surface(16) = {15};
Line Loop(17) = {3, 4, 1, 2};
Plane Surface(18) = {17};
Line Loop(19) = {12, -2, -11, -8};
Plane Surface(20) = {19};
Line Loop(21) = {7, 12, 3, -6};
Plane Surface(22) = {21};
Line Loop(23) = {9, -5, 1, -11};
Plane Surface(24) = {23};

Surface Loop(25) = {14, 22, 20, 18, 16, 24};

Volume(26) = {25};

Physical Volume("cube1") = {26};

//====================================
// Cube 2
//====================================
Point(9) = {-1+d, -1+d, -1+d, h};
Point(10) = {1-d, -1+d, -1+d, h};
Point(11) = {-1+d, 1-d, -1+d, h};
Point(12) = {1-d, 1-d, -1+d, h};
Point(13) = {1-d, 1-d, 1-d, h};
Point(14) = {1-d, -1+d, 1-d, h};
Point(15) = {-1+d, 1-d, 1-d, h};
Point(16) = {-1+d, -1+d, 1-d, h};

Line(13) = {11, 15};
Line(14) = {15, 13};
Line(15) = {13, 12};
Line(16) = {12, 11};
Line(17) = {11, 9};
Line(18) = {10, 12};
Line(19) = {10, 14};
Line(20) = {14, 16};
Line(21) = {16, 9};
Line(22) = {9, 10};
Line(23) = {16, 15};
Line(24) = {14, 13};

Line Loop(25) = {19, 20, 21, 22};
Plane Surface(26) = {25};
Line Loop(27) = {18, 16, 17, 22};
Plane Surface(28) = {27};
Line Loop(29) = {15, 16, 13, 14};
Plane Surface(30) = {29};
Line Loop(31) = {24, -14, -23, -20};
Plane Surface(32) = {31};
Line Loop(33) = {19, 24, 15, -18};
Plane Surface(34) = {33};
Line Loop(35) = {21, -17, 13, -23};
Plane Surface(36) = {35};

Surface Loop(37) = {26, 34, 32, 30, 28, 36, 14, 22, 20, 18, 16, 24};

Volume(38) = {37};

Physical Volume("cube2") = {38};

//====================================
// Cube 3
//====================================
Point(17) = {-1, -1, -1, h};
Point(18) = {1, -1, -1, h};
Point(19) = {-1, 1, -1, h};
Point(20) = {1, 1, -1, h};
Point(21) = {1, 1, 1, h};
Point(22) = {1, -1, 1, h};
Point(23) = {-1, 1, 1, h};
Point(24) = {-1, -1, 1, h};

Line(25) = {19, 23};
Line(26) = {23, 21};
Line(27) = {21, 20};
Line(28) = {20, 19};
Line(29) = {19, 17};
Line(30) = {18, 20};
Line(31) = {18, 22};
Line(32) = {22, 24};
Line(33) = {24, 17};
Line(34) = {17, 18};
Line(35) = {24, 23};
Line(36) = {22, 21};

Line Loop(37) = {31, 32, 33, 34};
Plane Surface(38) = {37};
Line Loop(39) = {30, 28, 29, 34};
Plane Surface(40) = {39};
Line Loop(41) = {27, 28, 25, 26};
Plane Surface(42) = {41};
Line Loop(43) = {36, -26, -35, -32};
Plane Surface(44) = {43};
Line Loop(45) = {31, 36, 27, -30};
Plane Surface(46) = {45};
Line Loop(47) = {33, -29, 25, -35};
Plane Surface(48) = {47};

Surface Loop(49) = {38, 46, 44, 42, 40, 48, 26, 34, 32, 30, 28, 36};

Volume(50) = {49};

Physical Volume("cube3") = {50};

Coherence;

