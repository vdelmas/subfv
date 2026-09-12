//+
Point(1) = {0, 0, 0, 1.0};
//+
Extrude {0, 1.5, 0} {
  Point{1}; 
}
//+
Extrude {0, 1.5, 0} {
  Point{2}; 
}
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Curve{2}; Curve{1}; 
}
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Curve{3}; Curve{7}; 
}
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Curve{10}; Curve{14}; 
}
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Curve{17}; Curve{21}; 
}
//+
Coherence;
//+
Transfinite Curve {12, 19, 26, 5, 4, 25, 18, 11} = 20 Using Progression 1;
//+
Transfinite Curve {3, 2, 17, 10} = 20 Using Progression 1;
//Transfinite Curve {3, 2, 17, 10} = 20 Using Progression 1.2;
//+
Transfinite Surface {13};
//+
Transfinite Surface {6};
//+
Transfinite Surface {27};
//+
Transfinite Surface {20};
//+
Recombine Surface {13, 6, 27, 20};
//+
Transfinite Curve {7, 1, 21, 14} = 20 Using Progression 1;
//Transfinite Curve {-7, -1, -21, -14} = 30 Using Progression 1.2;
//+
Recombine Surface {9, 30, 23, 16};
//+
nl=10;
Extrude {1, 0, 0} {
  Surface{6}; Surface{27}; Surface{30}; Surface{9}; Surface{20}; Surface{13}; Surface{16}; Surface{23}; Layers {4*nl}; Recombine;
}
//+
Extrude {6, 0, 0} {
  Surface{52}; Surface{74}; Surface{130}; Surface{152}; Surface{108}; Surface{91}; Surface{186}; Surface{169}; Layers {6*nl}; Recombine;
}
//+
Physical Volume("vol", 343) = {1, 2, 3, 4, 8, 7, 6, 5, 12, 11, 16, 15, 13, 14, 9, 10};
