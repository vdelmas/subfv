// Cartesian [0,400] x [-10,10] grid with dx=dy=2, extruded by one unit.
// The runner applies the tiny odd/even x perturbation to the y=0 node row.
SetFactory("Built-in");
Mesh.MshFileVersion = 4.1;
Mesh.ElementOrder = 1;
Mesh.RecombineAll = 1;
dx = 0.5;
Point(1) = {0, -10, 0, dx};
Point(2) = {400, -10, 0, dx};
Point(3) = {400, 10, 0, dx};
Point(4) = {0, 10, 0, dx};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Curve {1, 3} = 801;
Transfinite Curve {2, 4} = 41;
Transfinite Surface {1};
Recombine Surface {1};
out[] = Extrude {0, 0, 1} { Surface{1}; Layers{1}; Recombine; };
Physical Volume("fluid") = {out[1]};
Physical Surface("inflow") = {out[5]};
Physical Surface("outflow") = {out[3]};
Physical Surface("lower") = {out[2]};
Physical Surface("upper") = {out[4]};
