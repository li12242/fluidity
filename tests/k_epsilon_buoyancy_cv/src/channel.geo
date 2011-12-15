edgeLength = 0.01;

Point(1) = {0.0, 0.0, 0.0, 1};

Extrude {0.1, 0.0, 0.0} {
  Point{1}; Layers{0.1/edgeLength};
}

Extrude {0.0, 0.1, 0.0} {
  Line{1}; Layers{0.1/edgeLength};
}

Physical Line(1) = {1};  //reentrainment
Physical Line(2) = {2};  //lid
Physical Line(3) = {3};  //in
Physical Line(4) = {4};  //out
Physical Surface(6) = {5};
