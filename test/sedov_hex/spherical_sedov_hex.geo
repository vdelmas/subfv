nl = 31;
dl = 2.4;
Point(1) = {-1.2, -1.2, -1.2, 1.0};
Extrude {0, dl, 0} { Point{1}; Layers {nl}; Recombine; }
Extrude {dl, 0, 0} { Curve{1}; Layers {nl}; Recombine; }
Extrude {0, 0, dl} { Surface{5}; Layers {nl}; Recombine; }
Physical Volume('fluid', 28) = {1};
Physical Surface('side_surf', 29) = {14, 18, 22, 26};
