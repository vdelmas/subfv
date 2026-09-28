Lx = 4.0;
Ly = 2.0;
xc = 0.6;
yc = 1.0;
R = 0.25;

lc = 0.00463;

Point(1) = {0, 0, 0, lc};
Point(2) = {Lx, 0, 0, lc};
Point(3) = {Lx, Ly, 0, lc};
Point(4) = {0, Ly, 0, lc};

Point(5) = {xc, yc, 0, lc};
Point(6) = {xc+R, yc, 0, lc};
Point(7) = {xc, yc+R, 0, lc};
Point(8) = {xc-R, yc, 0, lc};
Point(9) = {xc, yc-R, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

Curve Loop(1) = {5, 6, 7, 8};
Curve Loop(2) = {1, 2, 3, 4};
Plane Surface(1) = {2, 1};

out1[] = Extrude {0, 0, 1e-3} { Surface{1}; Layers{1}; Recombine; };

Physical Volume("fluid", 1) = {out1[1]};
Physical Surface("bot_surf", 2) = Surface In BoundingBox{-1e-6, -1e-6, -1e-3-1e-6, Lx+1e-6, 1e-6, 1e-3+1e-6};
Physical Surface("top_surf", 3) = Surface In BoundingBox{-1e-6, Ly-1e-6, -1e-3-1e-6, Lx+1e-6, Ly+1e-6, 1e-3+1e-6};
Physical Surface("in_surf", 4) = Surface In BoundingBox{-1e-6, -1e-6, -1e-3-1e-6, 1e-6, Ly+1e-6, 1e-3+1e-6};
Physical Surface("out_surf", 5) = Surface In BoundingBox{Lx-1e-6, -1e-6, -1e-3-1e-6, Lx+1e-6, Ly+1e-6, 1e-3+1e-6};
Physical Surface("cyl_surf", 6) = Surface In BoundingBox{xc-R-1e-6, yc-R-1e-6, -1e-3-1e-6, xc+R+1e-6, yc+R+1e-6, 1e-3+1e-6};
