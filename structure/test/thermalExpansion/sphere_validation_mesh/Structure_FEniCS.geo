// Gmsh project created on Wed May 27 17:50:35 2020
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0, 0, 0, 0.19, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(2) = {0, 0, 0, 0.2, -Pi/2, Pi/2, 2*Pi};
//+
BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; }

//+
Box(3) = {-0.5, 0, -0.5, 1, -0.5, 1};
//+
BooleanDifference{ Volume{2}; Delete; }{ Volume{3}; Delete; }
//+
Box(3) = {-0.5, 0., 0, 1, 0.5, -0.5};
//+
BooleanDifference{ Volume{2}; Delete; }{ Volume{3}; Delete; }
//+
Box(3) = {0, 0., 0, -0.5, 0.5, 0.5};
//+
BooleanDifference{ Volume{2}; Delete; }{ Volume{3}; Delete; }
//+
Physical Surface("inner", 1) = {5};
//+
Physical Surface("outer", 2) = {1};
//+
Physical Surface("S1", 3) = {3};
//+
Physical Surface("S2", 4) = {4};
//+
Physical Surface("S3", 5) = {2};
//+
//Transfinite Curve {5, 6, 8} = 10 Using Progression 1;
//+
//Transfinite Curve {9, 2, 7, 1, 10, 3} = 20 Using Progression 1;
//+
//Transfinite Surface {4} = {3, 6, 4, 2};
//+
//Transfinite Surface {2} = {5, 1, 2, 4};
//+
//Transfinite Surface {3} = {5, 1, 3, 6};
//+
Physical Volume("volume", 6) = {2};
