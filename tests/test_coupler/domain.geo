Point(1) = {0, 0, 0, 1e+22};
Point(2) = {360, 0, 0, 1e+22};
Point(3) = {360, 360, 0, 1e+22};
Point(4) = {0, 360, 0, 1e+22};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {4, 1, 2, 3};
Plane Surface(6) = {6};
Physical Line(8) = {1, 2, 3, 4};
Physical Surface(7) = {6};
