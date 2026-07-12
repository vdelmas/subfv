insize=0.04;
outsize=0.08;
Point(1) = {-1, 0, 0, 1.0};
Extrude {-1.5, 0, 0} { Point{1}; }
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Curve{1}; }
Extrude {{0, 0, 1}, {0, 0, 0}, -Pi/2} { Curve{1}; }
Extrude {0, 0, 1e-3} { Surface{9}; Surface{5}; Layers {1}; Recombine; }
MeshSize {3, 1, 9, 21, 11, 32} = insize;
MeshSize {12, 2, 10, 17, 4, 28} = outsize;
Physical Surface("cyl_surf", 54) = {52, 30};
Physical Surface("in_surf", 55) = {44, 22};
Physical Surface("out_surf", 56) = {26, 48};
Physical Volume("fuild") = {1,2};
