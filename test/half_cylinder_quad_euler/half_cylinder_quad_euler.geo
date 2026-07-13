nlr=30;
nlt2=30;
Point(1) = {-1, 0, 0, 1.0};
//+
Extrude {-1.5, 0, 0} {
  Point{1}; Layers {nlr}; Recombine;
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Curve{1}; Layers{nlt2}; Recombine;
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, -Pi/2} {
  Curve{1}; Layers{nlt2}; Recombine;
}
//+
Extrude {0, 0, 1e-3} {
  Surface{9}; Surface{5}; Layers {1}; Recombine;
}
//+
Physical Surface("in_surf", 54) = {22, 44};
//+
Physical Surface("cyl_surf", 55) = {52, 30};
//+
Physical Surface("out_surf", 56) = {26, 48};
Physical Volume("fluid") = {1,2};

