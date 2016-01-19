// Gmsh project created on Sun Nov  8 22:29:17 2015
Point(1) = {10, 10, 0, 1.0};
Point(2) = {10, -10, 0, 1.0};
Point(3) = {-10, 10, 0, 1.0};
Point(4) = {-10, -10, 0, 1.0};
Point(5) = {0, 0, 0, 1.0};
Point(6) = {-2, 0, 0, 1.0};
Point(7) = {0, 2, 0, 1.0};
Point(8) = {2, 0, 0, 1.0};
Point(9) = {0, -2, 0, 1.0};
Line(1) = {3, 1};
Line(2) = {1, 2};
Line(3) = {2, 4};
Line(4) = {4, 3};
Circle(5) = {7, 5, 6};
Circle(6) = {6, 5, 9};
Circle(7) = {9, 8, 8};
Delete {
  Line{7};
}
Circle(7) = {7, 5, 8};
Circle(8) = {8, 5, 9};
Line Loop(9) = {4, 1, 2, 3};
Line Loop(10) = {5, 6, -8, -7};
Plane Surface(11) = {9, 10};
