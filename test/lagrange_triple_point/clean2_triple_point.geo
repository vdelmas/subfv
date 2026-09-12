//+
Point(1) = {0, 0, 0, 1.0};
//+
Extrude {0, 0.5, 0} {
  Point{1}; 
}
//+
Extrude {0, 1.0, 0} {
  Point{2}; 
}
//+
Extrude {0, 1.5, 0} {
  Point{3}; 
}
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Curve{3}; Curve{2}; Curve{1}; 
}
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Curve{4}; Curve{8}; Curve{12}; 
}
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Curve{15}; Curve{19}; Curve{23}; 
}
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Curve{26}; Curve{30}; Curve{34}; 
}
//+
Coherence;
//+
Transfinite Curve {6, 39, 28, 17, 5, 38, 27, 16, 9, 42, 31, 20} = 10 Using Progression 1; //theta
//+
Transfinite Curve {4, 3, 26, 15} = 10 Using Progression 1.0; //rext 1.2
//+
Transfinite Curve {-8, -2, -30, -19} = 10 Using Progression 1.0; //rin1 1.2
//+
Transfinite Surface {7};
//+
Transfinite Surface {40};
//+
Transfinite Surface {29};
//+
Transfinite Surface {18};
//+
Transfinite Surface {11};
//+
Transfinite Surface {44};
//+
Transfinite Surface {33};
//+
Transfinite Surface {22};
//+
Recombine Surface {14, 47, 36, 25};
//+
Transfinite Curve {1, 34, 23, 12} = 5 Using Progression 1; //rin2
//+
Recombine Surface {7, 40, 29, 18, 11, 44, 33, 22};
//+
nl=10;
Extrude {1, 0, 0} {
  Surface{7}; Surface{40}; Surface{29}; Surface{18}; Surface{22}; Surface{11}; Surface{44}; Surface{33}; Surface{36}; Surface{47}; Surface{14}; Surface{25}; Layers {4*nl}; Recombine;
}
//+
Extrude {6, 0, 0} {
  Surface{69}; Surface{91}; Surface{113}; Surface{135}; Surface{179}; Surface{201}; Surface{223}; Surface{157}; Surface{274}; Surface{257}; Surface{240}; Surface{291}; Layers {6*nl}; Recombine;
}
//+
Physical Volume("vol", 558) = {7, 6, 1, 14, 2, 3, 15, 8, 10, 18, 11, 23, 24, 12, 9, 21, 19, 5, 20, 22, 25, 17, 13, 4, 16};
