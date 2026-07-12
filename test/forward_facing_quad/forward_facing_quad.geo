// Forward Facing Step - structured quad mesh
// Domain: [0,xl]x[0,H] (inlet, full height) + [xl,L]x[h,H] (above step)
// Mach 3: rho=1.4, u=3, p=1, tmax=4.0

xl = 0.6; L = 3.0; h = 0.2; H = 1.0;
nl = 50;
dz = h / nl;

// 2D Points
Point(1) = {0,  0,  0};   // bottom-left
Point(2) = {0,  h,  0};   // inlet at step height
Point(3) = {0,  H,  0};   // top-left
Point(4) = {xl, 0,  0};   // step base
Point(5) = {xl, h,  0};   // step corner
Point(6) = {xl, H,  0};   // top at xl
Point(7) = {L,  h,  0};   // outlet bottom
Point(8) = {L,  H,  0};   // outlet top

// 2D Lines
Line(1)  = {1, 4};    // bottom floor  (y=0, x=[0,xl])      wall
Line(2)  = {4, 5};    // step face     (x=xl, y=[0,h])      wall
Line(3)  = {2, 5};    // interior      (y=h, x=[0,xl])      shared B/A
Line(4)  = {3, 6};    // top left      (y=H, x=[0,xl])      wall
Line(5)  = {1, 2};    // inlet lower   (x=0, y=[0,h])       freestream
Line(6)  = {2, 3};    // inlet upper   (x=0, y=[h,H])       freestream
Line(7)  = {5, 6};    // interior      (x=xl, y=[h,H])      shared A/C
Line(8)  = {6, 8};    // top right     (y=H, x=[xl,L])      wall
Line(9)  = {5, 7};    // step top      (y=h, x=[xl,L])      wall
Line(10) = {7, 8};    // outlet        (x=L, y=[h,H])       outflow

// 3 rectangular patches
// Patch B: [0,xl]x[0,h]  (inlet below step level)
Curve Loop(1) = {1, 2, -3, -5};
// Patch A: [0,xl]x[h,H]  (inlet above step level)
Curve Loop(2) = {3, 7, -4, -6};
// Patch C: [xl,L]x[h,H]  (above step, outlet region)
Curve Loop(3) = {9, 10, -8, -7};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

// Transfinite structured meshing
nx1 = 3*nl;   // x=[0,xl]  30 cells
nx2 = 12*nl;  // x=[xl,L] 120 cells
nz  = nl;     // y=[0,h]   10 cells (step height)
ny  = 4*nl;   // y=[h,H]   40 cells (above step)

Transfinite Curve{1, 3, 4}  = nx1+1;
Transfinite Curve{5, 2}     = nz+1;
Transfinite Curve{6, 7, 10} = ny+1;
Transfinite Curve{8, 9}     = nx2+1;

Transfinite Surface{1} = {1, 4, 5, 2};
Transfinite Surface{2} = {2, 5, 6, 3};
Transfinite Surface{3} = {5, 7, 8, 6};
Recombine Surface{1, 2, 3};

// Single-layer z extrusion
Extrude {0, 0, dz} {
  Surface{1}; Surface{2}; Surface{3};
  Layers{1}; Recombine;
}

// Physical groups: inlet (x=0) and outlet (x=L) only
// Wall surfaces default via is_wall(re==0)
// Surface IDs: 31=inlet_lower, 53=inlet_upper, 67=outlet, volumes 1+2+3
Physical Surface("inlet",  10) = {31, 53};
Physical Surface("outlet", 11) = {67};
Physical Volume("fluid",   12) = {1, 2, 3};
