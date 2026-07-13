// Annular mesh for 2D potential flow around a cylinder
// Inner cylinder radius r0 = 0.5, outer domain radius r1 = 5.5
// Triangular mesh in xy-plane, extruded 1 layer in z

r0 = 0.5;
r1 = 5.5;
ls_inner = 0.05;
ls_outer = 0.5;

// Center point
Point(1) = {0, 0, 0, 1.0};

// Points on inner circle (r0)
Point(2) = { r0,  0, 0, ls_inner};
Point(3) = {  0, r0, 0, ls_inner};
Point(4) = {-r0,  0, 0, ls_inner};
Point(5) = {  0,-r0, 0, ls_inner};

// Points on outer circle (r1)
Point(6) = { r1,  0, 0, ls_outer};
Point(7) = {  0, r1, 0, ls_outer};
Point(8) = {-r1,  0, 0, ls_outer};
Point(9) = {  0, -r1, 0, ls_outer};

// Inner circle arcs (clockwise = hole orientation)
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Outer circle arcs
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

// Annular surface: outer loop minus inner hole
Curve Loop(1) = {5, 6, 7, 8};
Curve Loop(2) = {1, 2, 3, 4};
Plane Surface(1) = {1, 2};

// Extrude one thin layer in z for the 2D simulation
Extrude {0, 0, 1.0} {
  Surface{1}; Recombine; Layers {1};
}

// Physical groups
// Only the outer boundary is tagged; inner cylinder uses default wall BC
Physical Surface("farfield", 100) = {21, 25, 29, 33};
Physical Volume("fluid", 101) = {1};
