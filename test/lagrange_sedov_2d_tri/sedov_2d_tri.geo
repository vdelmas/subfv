// Sedov 2D -- STRUCTURED triangle mesh (grid 24x24 over [-1.2,1.2]^2, each
// cell split into 2 triangles: transfinite, no Recombine on the surface).
// Same domain / init as lagrange_sedov_2d_quad, triangles instead of quads
// (as the Euler code's sedov_tri test).
// NOTE: the vertex-based Lagrangian scheme keeps a very small time step on
// triangles (collapsing central sub-cells), so this test uses a coarse mesh
// and a short t_max (0.15) -- just enough to compare the blast wave.
dl = 2.4;
Point(1) = {-1.2, -1.2, 0, 1.0};
MeshSize {1} = 0.02;
Extrude {0, dl, 0} { Point{1}; }
Extrude {dl, 0, 0} { Curve{1}; }
Extrude {0, 0, 1.0} { Surface{5}; Layers {1}; Recombine; }
Physical Volume("fluid", 28) = {1};
Mesh.Algorithm = 4;
