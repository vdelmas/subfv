// Annular mesh for 2D potential flow around a cylinder
// Inner cylinder radius r0 = 0.5, outer domain radius r1 = 10.0
// Structured quad mesh: 4 transfinite patches, extruded 1 layer in z

r0 = 0.5;
r1 = 10.0;
N_circ = 20;   // nodes per quarter arc (inner = outer, transfinite constraint)
N_rad  = 30;   // nodes along radial lines (packed near inner cylinder)

// Center point (needed for circle arcs)
Point(1) = {0, 0, 0, 1.0};

// Points on inner circle (r0)
Point(2) = { r0,  0, 0, 1.0};
Point(3) = {  0, r0, 0, 1.0};
Point(4) = {-r0,  0, 0, 1.0};
Point(5) = {  0,-r0, 0, 1.0};

// Points on outer circle (r1)
Point(6) = { r1,  0, 0, 1.0};
Point(7) = {  0, r1, 0, 1.0};
Point(8) = {-r1,  0, 0, 1.0};
Point(9) = {  0, -r1, 0, 1.0};

// Inner circle arcs
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Outer circle arcs
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

// Radial lines from inner to outer
Line(9)  = {2, 6};
Line(10) = {3, 7};
Line(11) = {4, 8};
Line(12) = {5, 9};

// Four patches of the annular domain
Curve Loop(1) = {1, 10, -5, -9};
Plane Surface(1) = {1};

Curve Loop(2) = {2, 11, -6, -10};
Plane Surface(2) = {2};

Curve Loop(3) = {3, 12, -7, -11};
Plane Surface(3) = {3};

Curve Loop(4) = {4, 9, -8, -12};
Plane Surface(4) = {4};

// Transfinite structured mesh
// Progression 1.1 packs nodes at the start of each radial line (inner cylinder)
Transfinite Curve {1, 2, 3, 4, 5, 6, 7, 8} = N_circ;
Transfinite Curve {9, 10, 11, 12} = N_rad Using Progression 1.1;

Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};

Recombine Surface {1, 2, 3, 4};

// Extrude one thin layer in z for 2D simulation
// out[i][0] = top surface, out[i][1] = volume, out[i][2..5] = lateral surfaces
// Lateral surface order follows the Curve Loop: {inner_arc, radial, outer_arc, radial}
// => out[i][4] is the outer (farfield) boundary surface for each patch
out1[] = Extrude {0, 0, 1e-3} { Surface{1}; Layers{1}; Recombine; };
out2[] = Extrude {0, 0, 1e-3} { Surface{2}; Layers{1}; Recombine; };
out3[] = Extrude {0, 0, 1e-3} { Surface{3}; Layers{1}; Recombine; };
out4[] = Extrude {0, 0, 1e-3} { Surface{4}; Layers{1}; Recombine; };

// Physical groups
// Only the outer boundary is tagged; inner cylinder uses default wall BC
Physical Surface("farfield", 100) = {out1[4], out2[4], out3[4], out4[4]};
Physical Volume("fluid", 101) = {out1[1], out2[1], out3[1], out4[1]};
