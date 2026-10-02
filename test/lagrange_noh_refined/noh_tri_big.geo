// Noh 2D, maillage triangulaire non structure sur le meme domaine [-1,1]^2 que
// noh_quad.geo, extrude d'une couche en z. lc choisi pour un nombre de mailles
// comparable au quad 100x100 : aire d'un triangle ~ 0.433*lc^2, donc pour 10000
// mailles sur une aire de 4, lc ~ 0.0304.
lc = 0.0304;

Point(1) = {-1.5, -1.5, 0, lc};
Point(2) = { 1.5, -1.5, 0, lc};
Point(3) = { 1.5,  1.5, 0, lc};
Point(4) = {-1.5,  1.5, 0, lc};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

// les indices des surfaces laterales sont recuperes de l'extrusion plutot que
// codes en dur : out[0] = face haute, out[1] = volume, out[2..5] = cotes
out[] = Extrude {0, 0, 1e-3} { Surface{1}; Layers{1}; Recombine; };

Physical Volume("fluid") = {out[1]};
Physical Surface("movingbound") = {out[2], out[3], out[4], out[5]};
Physical Surface("outer_wall") = {out[2], out[3], out[4], out[5]};
