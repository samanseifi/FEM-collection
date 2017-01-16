/*
** This is a simple finite element code for linear elastic material.
**
** 		Young's modulus: E
**		Poisson Ratio:   nu
**
*/

/* Currently makefile is problematic for matrix multiplication with Armadillo!
** For time being compile it with:
**
** g++ Main.cpp Mesh.cpp SolidElement.cpp -o fem1 -std=c++11 -larmadillo
**
*/

#include "Mesh.h"
#include "SolidElement.h"

#include <cmath>
#include <iostream>

#include <fstream>

#define pi 3.14159265358979323846

using namespace std;
using namespace arma;

vec shape(vec);
mat gradshape(vec);
mat set_B(mat);

int main(int argc, char *argv[]) {

  // Create a uniform mesh
  Mesh mesh(20, 20, 2, 2);

  // Plane-strain material tangent
  double E = 100.0;
  double nu = 0.3;

  // Initialize isotropic elastostatic material
  SolidElement Material(mesh, E, nu);
  Mat<double> K = Material.K;

  // Add nodal forces and boundary conditions
  vec f(2 * mesh.num_nodes(), fill::zeros);
  for (int i = 0; i < mesh.num_nodes(); i++) {
    if (mesh.nodes_set(i, 1) == 0.0) {
      K.row(2 * i) = zeros<mat>(1, K.n_cols);
      K.row(2 * i + 1) = zeros<mat>(1, K.n_cols);
      K(2 * i, 2 * i) = 1.0;
      K(2 * i + 1, 2 * i + 1) = 1.0;
    }
    if (mesh.nodes_set(i, 1) == mesh.ly) {
      double x = mesh.nodes_set(i, 0);
      double fbar = mesh.hx * (cos(8.0 * pi * x / mesh.lx)); // weird BC!
      fbar *= (x * (mesh.lx - x)) / ((mesh.lx) * (mesh.lx));
      f(2 * i + 1) = fbar;
      if (x == 0.0 || x == mesh.lx)
        f(2 * i + 1) *= 0.5;
    }
  }
  cout << "Solving linear system ..." << endl;
  mat U = solve(K, f);

  ofstream file1;
  ofstream file2;

  file1.open("u.dat");
  file1 << U << endl;
  file1.close();

  file2.open("mesh_info.dat");
  file2 << mesh.lx << " " << mesh.ly << " " << mesh.lx / mesh.hx + 1 << " "
        << mesh.ly / mesh.hy + 1 << endl;
  file2.close();

  return 0;
}
