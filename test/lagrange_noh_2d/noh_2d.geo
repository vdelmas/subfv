nl = 50;

Point(1) = {-1, -1, 0, 1.0};
Extrude {0, 2, 0} { Point{1}; Layers {nl}; Recombine; }
Extrude {2, 0, 0} { Curve{1}; Layers {nl}; Recombine; }
Extrude {0, 0, 1e-3} { Surface{5}; Layers {1}; Recombine; }

Physical Volume("fluid", 28) = {1};

// Outer walls: x=-1.2 (14), y=+1.2 (18), x=+1.2 (22), y=-1.2 (26)
Physical Surface("outer_wall", 1) = {14, 18, 22, 26};
//+
Physical Surface("movingbound", 29) = {18, 14, 22, 26};
