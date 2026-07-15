//+
Point(1) = {-0.5, -0.5, -0.5, 1.0};
//+
Extrude {1, 0, 0} {
  Point{1}; Layers {50}; Recombine;
}
//+
Extrude {0, 1, 0} {
  Curve{1}; Layers {50}; Recombine;
}
//+
Extrude {0, 0, 1} {
  Surface{5}; Layers {1}; Recombine;
}
//+
Physical Volume("fluid", 28) = {1};
//+
Physical Surface("all_surf", 29) = {26, 5, 22, 18, 27, 14};
