Point(1) = {0.0, 0.0, 0.0, 300};
Point(2) = {0.0, 10000.0, 0.0, 300};
Point(3) = {40000.0, 10000.0, 0.0, 300};
Point(4) = {40000.0, 0.0, 0.0, 300};
Point(5) = {20500.0, 0.0, 0.0, 300};
Point(6) = {20000.0, 0.0, 0.0, 300};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line Loop(7) = {1,2,3,4,5,6};

Plane Surface(8) = {7};

// Left side
Physical Line(111) = {1};
// Right side
Physical Line(222) = {3};
// Top
Physical Line(333) = {2};
// Bottom (EXCLUDING VENT)
Physical Line(444) = {4, 6};
// Bottom (VENT)
Physical Line(999) = {5};

// This is just to ensure all the interior
// elements get written out. 
Physical Surface(9) = {8};
