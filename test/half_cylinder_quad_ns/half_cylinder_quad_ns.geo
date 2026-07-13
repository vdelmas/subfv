nlr=64;//128
nlthet=16;
Point(1) = {0, -1., 0, 1.0};
Point(2) = {0, -1.7, 0, 1.0};
Line(1) = {1, 2};
Transfinite Curve {1} = nlr Using Beta_HWall 2e-5;
Extrude {{0, 0, 1}, {0, 0, 0}, -Pi/4} { Curve{1}; Layers{nlthet}; Recombine; }
Extrude {{0, 0, 1}, {0, 0, 0}, -Pi/4} { Curve{2}; Layers{nlthet}; Recombine; }
Extrude {{0, 0, 1}, {0, 0, 0}, -Pi/4} { Curve{6}; Layers{nlthet}; Recombine; }
Extrude {{0, 0, 1}, {0, 0, 0}, -Pi/4} { Curve{10}; Layers{nlthet}; Recombine; }
Extrude {0, 0, 1e-3} { Surface{9}; Surface{5}; Surface{17}; Surface{13}; Layers {1}; Recombine; }
Physical Surface('surf_in') = {96, 30, 52, 74};
Physical Surface('surf_cyl') = {104, 82, 38, 60};
Physical Surface('surf_out') = {48, 78};
Physical Volume('fluid') = {1, 2, 3, 4};
