// Square [-L/2, L/2]^2 extruded by one layer of thickness dz.
// Mesh size and element type are set from the command line:
//   gmsh -3 vortex.geo -setnumber ls 0.25 -setnumber quad 0 -o mesh.msh
DefineConstant[
  ls   = {0.3, Name "ls"},     // target edge length
  L    = {10.0, Name "L"},    // domain side: [-5,5] x [-5,5]
  dz   = {0.1, Name "dz"},     // extrusion thickness
  quad = {0, Name "quad"}      // 0: triangles, 1: quads (transfinite)
];

Point(1) = {-L/2, -L/2, 0, ls};
Point(2) = { L/2, -L/2, 0, ls};
Point(3) = { L/2,  L/2, 0, ls};
Point(4) = {-L/2,  L/2, 0, ls};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

If (quad == 1)
  n = Round(L/ls) + 1;
  Transfinite Curve{1, 2, 3, 4} = n;
  Transfinite Surface{1};
  Recombine Surface{1};
EndIf

out[] = Extrude {0, 0, dz} { Surface{1}; Layers{1}; Recombine; };
// out[0] = top, out[1] = volume, out[2..5] = sides built on lines 1..4
Physical Volume("fluid") = {out[1]};
Physical Surface("far_surf") = {out[2], out[3], out[4], out[5]};
