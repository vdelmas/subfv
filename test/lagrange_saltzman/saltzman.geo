nx = 100; ny = 10;

Point(1) = {0,   0,   0, 1.0};
Point(2) = {1,   0,   0, 1.0};
Point(3) = {1,   0.1, 0, 1.0};
Point(4) = {0,   0.1, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Curve{1, 3} = nx + 1;
Transfinite Curve{2, 4} = ny + 1;
Transfinite Surface{1};
Recombine Surface{1};

out[] = Extrude {0, 0, 0.1} { Surface{1}; Layers{1}; Recombine; };
// out[0] = top surface (z=0.1)
// out[1] = volume
// out[2] = side from Curve{1} (y=0)
// out[3] = side from Curve{2} (x=1, right wall)
// out[4] = side from Curve{3} (y=0.1)
// out[5] = side from Curve{4} (x=0, left wall = piston)

Physical Volume("fluid", 28)       = {out[1]};
Physical Surface("right_wall", 29) = {out[3]};
Physical Surface("left_wall", 30)  = {out[5]};
