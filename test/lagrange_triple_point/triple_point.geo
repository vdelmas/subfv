// Maillage 4x plus fin dans x=[0,1] (zone haute pression) que dans x=[1,7]
// Resolution: h_x ~ 1/28 dans [0,1], h_x ~ 6/42 dans [1,7]
nx1 = 28;  // cellules dans x=[0,1]
nx2 = 42;  // cellules dans x=[1,7]
ny  = 48;  // cellules en y=[-3,3]
nz  = 48;  // cellules en z=[-3,3]

// Points dans le plan z=-3
Point(1) = {0,  -3, -3, 1.0};
Point(2) = {1,  -3, -3, 1.0};
Point(3) = {7,  -3, -3, 1.0};
Point(4) = {7,   3, -3, 1.0};
Point(5) = {1,   3, -3, 1.0};
Point(6) = {0,   3, -3, 1.0};

// Courbes bloc gauche [0,1]x[-3,3]
Line(1) = {1, 2};   // bas   x=[0,1]
Line(2) = {2, 5};   // droit x=1
Line(3) = {5, 6};   // haut  x=[0,1]
Line(4) = {6, 1};   // gauche x=0

// Courbes bloc droit [1,7]x[-3,3]
Line(5) = {2, 3};   // bas   x=[1,7]
Line(6) = {3, 4};   // droit x=7
Line(7) = {4, 5};   // haut  x=[1,7]
// -Line(2) = interface partagee x=1 (utilisee en sens inverse)

Curve Loop(1) = {1,  2,  3, 4};
Curve Loop(2) = {5,  6,  7, -2};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Transfinite Curve{1,  3} = nx1 + 1;
Transfinite Curve{4,  2} = ny  + 1;
Transfinite Surface{1} = {1, 2, 5, 6};
Recombine Surface{1};

Transfinite Curve{5,  7} = nx2 + 1;
Transfinite Curve{6,  2} = ny  + 1;
Transfinite Surface{2} = {2, 3, 4, 5};
Recombine Surface{2};

out1[] = Extrude {0, 0, 6} { Surface{1}; Layers{nz}; Recombine; };
out2[] = Extrude {0, 0, 6} { Surface{2}; Layers{nz}; Recombine; };

Physical Volume("fluid", 28) = {out1[1], out2[1]};
