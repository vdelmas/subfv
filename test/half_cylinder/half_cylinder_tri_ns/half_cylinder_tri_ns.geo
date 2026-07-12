r1 = 1.04;
r2 = 1.7;
nlr=16;
nlthet=32;
ls2 = 3.1415*r1/(4*nlthet);
ls3 = 3.1415*r2/(4*nlthet);
Point(1) = {0, -1., 0, 1.0};
Point(2) = {0, -r1, 0, ls2};
Point(3) = {0, -r2, 0, ls3};
Line(1) = {1, 2};
Line(2) = {2, 3};
Transfinite Curve {1} = nlr Using Beta_HWall 1e-4;
Extrude {{0, 0, 1}, {0, 0, 0},-Pi/4} { Curve{1}; Layers{nlthet}; Recombine; }
Extrude {{0, 0, 1}, {0, 0, 0},-Pi/4} { Curve{3}; Layers{nlthet}; Recombine; }
Extrude {{0, 0, 1}, {0, 0, 0},-Pi/4} { Curve{7}; Layers{nlthet}; Recombine; }
Extrude {{0, 0, 1}, {0, 0, 0},-Pi/4} { Curve{11}; Layers{nlthet}; Recombine; }
Extrude {{0, 0, 1}, {0, 0, 0},-Pi/4} { Curve{2}; }
Extrude {{0, 0, 1}, {0, 0, 0},-Pi/4} { Curve{19}; }
Extrude {{0, 0, 1}, {0, 0, 0},-Pi/4} { Curve{23}; }
Extrude {{0, 0, 1}, {0, 0, 0},-Pi/4} { Curve{27}; }
Coherence;
//+
Extrude {0, 0, 1e-3} {
  Surface{30}; Surface{26}; Surface{34}; Surface{18}; Surface{14}; Surface{10}; Surface{6}; Surface{22}; Layers {1}; Recombine;
}
//+
Physical Surface("surf_in", 211) = {91, 47, 69, 201};
//+
Physical Surface("surf_cyl", 212) = {121, 143, 165, 187};
//+
Physical Surface("surf_out", 213) = {197, 95, 175, 117};
//+
Physical Volume("fluid", 214) = {8, 7, 6, 2, 5, 1, 4, 3};
