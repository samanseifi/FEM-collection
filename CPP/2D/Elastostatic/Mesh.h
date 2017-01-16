#ifndef MESH_H_
#define MESH_H_

#include <armadillo>

using namespace arma;

class Mesh {

public:
  Mesh(int, int, double, double);
  virtual ~Mesh();
  int num_nodes();
  int num_elements();

  mat nodes_set; // set of nodes coordinates
  mat conns_set; // set of connectivities by node numbers

  double lx, ly; // length in x and y direction
  double hx, hy; // x and y spacings

private:
  int nx;
  int ny;
};

#endif
