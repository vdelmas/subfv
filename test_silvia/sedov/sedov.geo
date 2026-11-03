// 2D Cartesian Sedov domain [-1.2,1.2]^2, extruded by one unit in z.
// 121 cells per side leave one unique cell centered at the origin.
SetFactory("Built-in");
Mesh.MshFileVersion = 4.1;
Mesh.ElementOrder = 1;
Mesh.RecombineAll = 1;
n = 121;
dx = 2.4/n;

Point(1) = {-1.2, -1.2, 0, dx};
Point(2) = { 1.2, -1.2, 0, dx};
Point(3) = { 1.2,  1.2, 0, dx};
Point(4) = {-1.2,  1.2, 0, dx};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Curve {1, 3} = n+1;
Transfinite Curve {2, 4} = n+1;
Transfinite Surface {1};
Recombine Surface {1};

out[] = Extrude {0, 0, 1} { Surface{1}; Layers{1}; Recombine; };
Physical Volume("fluid") = {out[1]};
Physical Surface("walls") = {out[2], out[3], out[4], out[5]};
