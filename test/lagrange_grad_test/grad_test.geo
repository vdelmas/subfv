// Unit square [0,1]^2 extruded dz=1e-3, unstructured triangle base (prism 3D)
// Gradient test: p = 1 + x + y, v = 0, rho = 1
// Exact grad_p = [1, 1, 0], exact grad_v = 0, exact div_v = 0
lc = 0.01;
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
// No Recombine on Surface -> triangle base -> prism 3D elements
ext[] = Extrude {0, 0, 1e-3} {
  Surface{1}; Layers{1}; Recombine;
};
Physical Volume("fluid", 1) = {ext[1]};
Physical Surface("wall", 2) = {1, ext[0], ext[2], ext[3], ext[4], ext[5]};
