n = 1;
Point(1) = {-1, -1, 0, 1.0};
Extrude {2, 0, 0} { Point{1}; Layers {n*200}; Recombine; }
Extrude {0, 2, 0} { Curve{1}; Layers {200}; Recombine; }
Extrude {0, 0, 1e-3} { Surface{5}; Layers {1}; Recombine; }
Physical Volume("fluid", 28) = {1};
