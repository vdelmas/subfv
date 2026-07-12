ls = 0.02;
Point(1) = {0, 0, 0, 1.0};
//+
Extrude {1, 0, 0} {
  Point{1};
}
//+
Extrude {0, 1, 0} {
  Curve{1};
}
//+
Extrude {0, 0, 1} {
  Surface{5}; Layers {1}; Recombine;
}
//+
Physical Volume("fluid", 28) = {1};
//+
MeshSize {3, 4, 10, 14, 1, 5, 2, 6} = ls;
