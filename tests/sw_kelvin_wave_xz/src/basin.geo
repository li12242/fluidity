r=250000.0;
Point(1) = {0,0,0,100000.0};
Point(2) = {r,0,0,10000.0};
Point(3) = {0,0,r,10000.0};
Point(4) = {-r,0,0,10000.0};
Point(5) = {0,0,-r,10000.0};
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Physical Line(3) = {1,2,3,4};
Line Loop(4) = {1,2,3,4};
Plane Surface(5) = {4};
Physical Surface(6) = {5};
