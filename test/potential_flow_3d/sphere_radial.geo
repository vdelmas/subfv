// Spherical shell mesh for 3D potential flow around a sphere
// Inner sphere radius r = 0.5, outer domain radius r2 = 5.5
// Structured hex mesh (radial extrusion from inner sphere)
// Only sphere_out is tagged; sphere_in uses default wall BC.

nl = 10;
ls = 2./nl;

r = 0.5;
r2 = 5.5;

Point(1000) = {0, 0, 0, ls};

Point(1) = {r, 0, 0, ls};
Point(2) = {0, r, 0, ls};
Point(3) = {-r, 0, 0, ls};
Point(4) = {0, -r, 0, ls};
Point(5) = {0, 0, r, ls};
Point(6) = {0, 0, -r, ls};
Circle(1) = {3, 1000, 5};
Circle(2) = {5, 1000, 1};
Circle(3) = {1, 1000, 6};
Circle(4) = {6, 1000, 3};
Circle(5) = {3, 1000, 2};
Circle(6) = {2, 1000, 1};
Circle(7) = {2, 1000, 5};
Circle(8) = {5, 1000, 4};
Circle(9) = {4, 1000, 3};
Circle(10) = {6, 1000, 4};
Circle(11) = {4, 1000, 1};
Circle(12) = {6, 1000, 2};

Curve Loop(1) = {7, -1, 5};
Surface(1) = {1};
Curve Loop(2) = {7, 2, -6};
Surface(2) = {2};
Curve Loop(3) = {6, 3, 12};
Surface(3) = {3};
Curve Loop(4) = {12, -5, -4};
Surface(4) = {4};
Curve Loop(5) = {3, 10, 11};
Surface(5) = {5};
Curve Loop(6) = {10, 9, -4};
Surface(6) = {6};
Curve Loop(7) = {9, 1, 8};
Surface(7) = {7};
Curve Loop(8) = {11, -2, 8};
Surface(8) = {8};
Extrude{ Surface{-1,2,3,-4,-5,6,-7,8}; Layers{nl,r2-r}; Recombine; }

Physical Surface('sphere_out', 151) = {114, 148, 97, 80, 46, 29, 131, 63};
Physical Volume('fluid', 152) = {1, 2, 8, 7, 6, 5, 3, 4};
