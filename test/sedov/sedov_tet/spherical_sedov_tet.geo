nl = 31;
ls = 2.4/nl;
dl = 2.4;
MeshSize {1} = ls;
Point(1) = {-1.2, -1.2, -1.2, 1.0};
Extrude {0, dl, 0} { Point{1};}
Extrude {dl, 0, 0} { Curve{1};}
Extrude {0, 0, dl} { Surface{5};}
Physical Volume('fluid', 28) = {1};
Physical Surface('side_surf', 29) = {14, 18, 22, 26};
MeshSize {2, 4, 10, 6, 1, 3, 14, 5} = ls;
