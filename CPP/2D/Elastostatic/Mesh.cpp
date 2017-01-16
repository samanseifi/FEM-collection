#include <armadillo>
#include <iostream>

#include "Mesh.h"

using namespace std;
using namespace arma;

Mesh::~Mesh() {}

/*
** Generating a uniform mesh
*
** (Note: Initilizing Armadillo's matrices should be done directly)
*/
Mesh::Mesh(int ex, int ey, double l_x, double l_y)
    : lx(l_x), ly(l_y), hx(l_x / ex), hy(l_y / ey),
      nodes_set((ex + 1) * (ey + 1), 2, fill::zeros),
      conns_set(ex * ey, 4, fill::zeros) {

  cout << "Constructing Mesh..." << endl;

  nx = ex + 1; // number of nodes in x
  ny = ey + 1; // number of nodes in y

  // Discretizing the lengths lx and ly
  vec x = linspace<vec>(0, lx, nx);
  vec y = linspace<vec>(0, ly, ny);

  // Creating the nodes set coordinates (x,y)
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      nodes_set(i + j * nx, 0) = x(i);
      nodes_set(i + j * nx, 1) = y(j);
    }
  }

  // Creating the connectivity matrix
  for (int j = 0; j < ey; j++) {
    for (int i = 0; i < ex; i++) {
      vec ele_conn(4, fill::zeros);

      // Numbering the nodes of the element i-th
      ele_conn(0) = i + j * nx;          // NORTH WEST
      ele_conn(1) = i + j * nx + 1;      // NORTH EAST
      ele_conn(2) = i + j * nx + 1 + nx; // SOUTH EAST
      ele_conn(3) = i + j * nx + nx;     // SOUTH WEST

      // Store to the global connectivity matrix
      conns_set.row(i + j * ex) = floor(trans(ele_conn));
    }
  }
}

/*
** Return the number of nodes
*/
int Mesh::num_nodes() { return nx * ny; }

/*
** Return the number of elements
*/
int Mesh::num_elements() { return (nx - 1) * (ny - 1); }
