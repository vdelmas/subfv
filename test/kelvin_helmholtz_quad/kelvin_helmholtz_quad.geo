// Kelvin-Helmholtz instability — structured quad mesh
// Domain: [-0.5, 2.5] x [-0.5, 0.5]

lx   = 3.0;
nl_y = 50;
nl_x = 3*nl_y;
dz   = 1.0 / nl_y;

Point(1) = {-0.5, -0.5, 0, 1.0};
Extrude {0, 1., 0}  { Point{1}; Layers {nl_y}; Recombine; }
Extrude {lx, 0, 0}  { Curve{1}; Layers {nl_x}; Recombine; }
Extrude {0, 0, dz}  { Surface{5}; Layers {1}; Recombine; }

Physical Volume("fluid",  28) = {1};
Physical Surface("left",  29) = {14};
Physical Surface("right", 30) = {22};
Physical Surface("topbot",31) = {18, 26};
