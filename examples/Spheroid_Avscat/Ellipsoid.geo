//
// gmsh geometry specification for a sphere of radius R=1
// 

//************************************************************
//* input parameters      
//************************************************************
R = 1.0;    // main radius
// Rsphere: radius of surface-equivalent sphere

// --------------- (D/L = 0.5, Rsphere = 1.0) --------------
//ax = 1.2032*R;   // x-axis
//ay = 1.2032*R;   // y-axis
//az = 0.6016*R;   // z-axis

// --------------- (D/L = 2.0, Rsphere = 1.0) --------------
ax = 0.7652*R;   // x-axis
ay = 0.7652*R;   // y-axis
az = 1.5304*R;   // z-axis

// center of ellipsoid
cx = 0.0;
cy = 0.0;
cz = 0.0;

//************************************************************
//* meshing finenesses ***************************************
//************************************************************
lx = 0.30*ax;  // fineness at north pole
ly = 0.30*ay;  // fineness at equator
lz = 0.15*az;  // fineness at south pole

////////////////////////////////////////
// ellipsoid

// middle point
Point(9) = {cx, cy, cz};

// major axis
Point(10) = {  ax + cx, 0.0 + cy, 0.0 + cz, lx};
Point(11) = { 0.0 + cx, 0.0 + cy,  az + cz, lz};
Point(12) = { -ax + cx, 0.0 + cy, 0.0 + cz, lx};
Point(13) = { 0.0 + cx, 0.0 + cy, -az + cz, lz};
Point(14) = { 0.0 + cx,  ay + cy, 0.0 + cz, ly};
Point(15) = { 0.0 + cx, -ay + cy, 0.0 + cz, ly};

// ellipsoid, lines
Ellipse(13) = {10, 9, 11, 11};
Ellipse(14) = {12, 9, 11, 11};
Ellipse(15) = {12, 9, 13, 13};
Ellipse(16) = {13, 9, 10, 10};

Ellipse(17) = {14, 9, 10, 10};
Ellipse(18) = {14, 9, 11, 11};
Ellipse(19) = {14, 9, 12, 12};
Ellipse(20) = {14, 9, 13, 13};

Ellipse(21) = {15, 9, 10, 10};
Ellipse(22) = {15, 9, 11, 11};
Ellipse(23) = {15, 9, 12, 12};
Ellipse(24) = {15, 9, 13, 13};


// line loops
Line Loop(1) = {16,-17,20};
Ruled Surface(101) = {1};
Line Loop(2) = {20,-15,-19};
Ruled Surface(102) = {2};
Line Loop(3) = {15,-24,23};
Ruled Surface(103) = {3};
Line Loop(4) = {16,-21,24};
Ruled Surface(104) = {4};
Line Loop(5) = {23,14,-22};
Ruled Surface(105) = {5};
Line Loop(6) = {21,13,-22};
Ruled Surface(106) = {6};
Line Loop(7) = {14,-18,19};
Ruled Surface(107) = {7};
Line Loop(8) = {18,-13,-17};
Ruled Surface(108) = {8};
Physical Surface(1) = {101,102,103,104,105,106,107,108};

//************************************************************
//* reference point to get outward-pointing surface normals right
//************************************************************
Physical Point(1) = {9};
