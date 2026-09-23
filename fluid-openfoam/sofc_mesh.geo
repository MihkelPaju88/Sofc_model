// Gmsh script for 2D Axisymmetric Microtubular SOFC (5-degree wedge)

// Parameters (in meters)
L = 0.05;      // Length
r_fuel = 0.001;  // Fuel channel radius (1mm)
r_porous = 0.0015; // Porous layer radius (1.5mm)
r_air = 0.003;   // Air channel radius (3mm)

// Grid densities
n_x = 50;      // Cells along length
n_r_fuel = 20; // Cells in fuel channel
n_r_porous = 10; // Cells in porous layer
n_r_air = 20;  // Cells in air channel

// Angle for OpenFOAM wedge (2.5 degrees on each side of Y-axis = 5 deg total)
theta = 2.5 * Pi / 180;

// Points (X, Y, Z) - Front Face
Point(1) = {0, 0, 0};
Point(2) = {L, 0, 0};

Point(3) = {0, r_fuel * Cos(theta), -r_fuel * Sin(theta)};
Point(4) = {L, r_fuel * Cos(theta), -r_fuel * Sin(theta)};

Point(5) = {0, r_porous * Cos(theta), -r_porous * Sin(theta)};
Point(6) = {L, r_porous * Cos(theta), -r_porous * Sin(theta)};

Point(7) = {0, r_air * Cos(theta), -r_air * Sin(theta)};
Point(8) = {L, r_air * Cos(theta), -r_air * Sin(theta)};

// Points (X, Y, Z) - Back Face
Point(9) = {0, r_fuel * Cos(theta), r_fuel * Sin(theta)};
Point(10) = {L, r_fuel * Cos(theta), r_fuel * Sin(theta)};

Point(11) = {0, r_porous * Cos(theta), r_porous * Sin(theta)};
Point(12) = {L, r_porous * Cos(theta), r_porous * Sin(theta)};

Point(13) = {0, r_air * Cos(theta), r_air * Sin(theta)};
Point(14) = {L, r_air * Cos(theta), r_air * Sin(theta)};

// Lines
Line(1) = {1, 2}; // Axis

// Fuel Channel Front
Line(2) = {3, 4};
Line(3) = {1, 3};
Line(4) = {2, 4};

// Porous Layer Front
Line(5) = {5, 6};
Line(6) = {3, 5};
Line(7) = {4, 6};

// Air Channel Front
Line(8) = {7, 8};
Line(9) = {5, 7};
Line(10) = {6, 8};

// Fuel Channel Back
Line(11) = {9, 10};
Line(12) = {1, 9};
Line(13) = {2, 10};

// Porous Layer Back
Line(14) = {11, 12};
Line(15) = {9, 11};
Line(16) = {10, 12};

// Air Channel Back
Line(17) = {13, 14};
Line(18) = {11, 13};
Line(19) = {12, 14};

// Connect Front to Back (Arcs for perfectly cylindrical geometry)
Circle(20) = {3, 1, 9}; // Inlet Fuel
Circle(21) = {4, 2, 10}; // Outlet Fuel
Circle(22) = {5, 1, 11}; // Inlet Porous
Circle(23) = {6, 2, 12}; // Outlet Porous
Circle(24) = {7, 1, 13}; // Inlet Air
Circle(25) = {8, 2, 14}; // Outlet Air

// Transfinite Lines (Setting grid densities)
Transfinite Line {1, 2, 5, 8, 11, 14, 17} = n_x + 1; // Length
Transfinite Line {3, 4, 12, 13} = n_r_fuel + 1; // Radial Fuel
Transfinite Line {6, 7, 15, 16} = n_r_porous + 1; // Radial Porous
Transfinite Line {9, 10, 18, 19} = n_r_air + 1; // Radial Air
Transfinite Line {20, 21, 22, 23, 24, 25} = 2; // 1 cell thick in theta

// Surfaces - Front
Line Loop(1) = {1, 4, -2, -3}; Surface(1) = {1}; // Fuel
Line Loop(2) = {2, 7, -5, -6}; Surface(2) = {2}; // Porous
Line Loop(3) = {5, 10, -8, -9}; Surface(3) = {3}; // Air

// Surfaces - Back
Line Loop(4) = {1, 13, -11, -12}; Surface(4) = {4}; // Fuel
Line Loop(5) = {11, 16, -14, -15}; Surface(5) = {5}; // Porous
Line Loop(6) = {14, 19, -17, -18}; Surface(6) = {6}; // Air

// Surfaces - Inlet
Line Loop(7) = {3, 20, -12}; Surface(7) = {7}; // Fuel
Line Loop(8) = {6, 22, -15, -20}; Surface(8) = {8}; // Porous
Line Loop(9) = {9, 24, -18, -22}; Surface(9) = {9}; // Air

// Surfaces - Outlet
Line Loop(10) = {4, 21, -13}; Surface(10) = {10}; // Fuel
Line Loop(11) = {7, 23, -16, -21}; Surface(11) = {11}; // Porous
Line Loop(12) = {10, 25, -19, -23}; Surface(12) = {12}; // Air

// Surfaces - Interfaces
Line Loop(13) = {2, 21, -11, -20}; Surface(13) = {13}; // Fuel/Porous
Line Loop(14) = {5, 23, -14, -22}; Surface(14) = {14}; // Porous/Air

// Surface - Outer Wall
Line Loop(15) = {8, 25, -17, -24}; Surface(15) = {15};

// Transfinite Surfaces
Transfinite Surface {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
Recombine Surface {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

// Volumes
Surface Loop(1) = {1, 4, 7, 10, 13}; Volume(1) = {1}; // Fuel
Surface Loop(2) = {2, 5, 8, 11, 13, 14}; Volume(2) = {2}; // Porous (Note: 13 shared)
Surface Loop(3) = {3, 6, 9, 12, 14, 15}; Volume(3) = {3}; // Air (Note: 14 shared)

Transfinite Volume {1, 2, 3};

// Physical Groups (For OpenFOAM mapping)
Physical Volume("fuel_channel") = {1};
Physical Volume("porous_layer") = {2};
Physical Volume("air_channel") = {3};

Physical Surface("inlet") = {7, 8, 9};
Physical Surface("outlet") = {10, 11, 12};
Physical Surface("axis") = {}; // OpenFOAM handles this automatically if Y=0
Physical Surface("front") = {1, 2, 3};
Physical Surface("back") = {4, 5, 6};
Physical Surface("outer_wall") = {15};
