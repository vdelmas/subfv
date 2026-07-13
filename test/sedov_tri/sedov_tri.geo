nl = 121;
dl = 2.4;
ls = 2.4/nl;
Point(1) = {-1.2, -1.2, -1.2, 1.0};
MeshSize {1} = ls;
Extrude {0, dl, 0} { Point{1}; }
Extrude {dl, 0, 0} { Curve{1}; }
Extrude {0, 0, 1.0} { Surface{5}; Layers {1}; Recombine; }
Physical Volume('fluid', 28) = {1};
Physical Surface('side_surf', 29) = {14, 18, 22, 26};
