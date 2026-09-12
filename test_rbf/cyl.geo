lsin=0.2;
lsout=20*lsin;
//+
Point(1) = {-1, 0, 0, 1.0};
//+
Extrude {-10, 0, 0} {
  Point{1}; 
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Curve{1}; 
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Curve{2}; 
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Curve{6}; 
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Curve{10}; 
}
//+
Extrude {0, 0, 1} {
  Surface{17}; Surface{13}; Surface{9}; Surface{5}; Layers {1}; Recombine;
}
//+
Physical Surface("in_cyl", 106) = {104, 38, 60, 82};
//+
Physical Volume("fluid", 107) = {1, 2, 3, 4};
//+
MeshSize {1, 23, 13, 11, 9, 24, 26, 3} = lsin;
//+
MeshSize {19, 2, 14, 12, 25, 10, 27, 4} = lsout;
