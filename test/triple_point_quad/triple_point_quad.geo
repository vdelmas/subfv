// Triple-point shock interaction -- structured quad mesh
// Domain: [0,7] x [0,3], single material (gamma = 1.4)
// Reference resolution: 700 x 300 (square cells, h = 0.01)

lx   = 7.0;
ly   = 3.0;
nl_x = 700;
nl_y = 300;
dz   = ly / nl_y;

Point(1) = {0, 0, 0, 1.0};
Extrude {0, ly, 0}  { Point{1}; Layers {nl_y}; Recombine; }
Extrude {lx, 0, 0}  { Curve{1}; Layers {nl_x}; Recombine; }
Extrude {0, 0, dz}  { Surface{5}; Layers {1}; Recombine; }

Physical Volume("fluid",  28) = {1};
Physical Surface("walls",  29) = {14, 18, 22, 26};
