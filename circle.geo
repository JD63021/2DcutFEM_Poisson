// Gmsh project created on Fri Sep 19 12:39:51 2025
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 10, 0, 2*Pi};
//+
Physical Curve(2) = {1};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
Physical Surface(3) = {1};
