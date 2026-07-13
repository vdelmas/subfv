//+
Point(1) = {0, 0, 0, 1.0};
//+
Extrude {0, 0, 0.1} {
  Point{1}; Layers {10}; Recombine;
}
//+
Extrude {0, 0.1, 0} {
  Curve{1}; Layers {10}; Recombine;
}
//+
Extrude {0.5, 0, 0} {
  Surface{5}; Layers {50}; Recombine;
}
//+
Physical Volume("fluid", 28) = {1};
//+
Physical Surface("left_surf", 29) = {5};
//+
Physical Surface("right_surf", 30) = {27};
