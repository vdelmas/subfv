nl=61;
Point(1) = {-1.2, -1.2, 0, 1.0};
Extrude {0, 2.4, 0} { Point{1}; Layers {nl}; Recombine; }
Extrude {2.4, 0, 0} { Curve{1}; Layers {nl}; Recombine; }
Extrude {0, 0, 1.0} { Surface{5}; Layers {1}; Recombine; }
Physical Volume("fluid", 28) = {1};
