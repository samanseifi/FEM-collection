#include "SolidElement.h"

#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

SolidElement::~SolidElement() {}

SolidElement::SolidElement(Mesh mesh, double E, double nu)
    : K(2 * mesh.num_nodes(), 2 * mesh.num_nodes(), fill::zeros) {

  Mat<double> D(3, 3, fill::zeros);
  D(0, 0) = D(1, 1) = 1.0 - nu;
  D(0, 1) = D(1, 0) = nu;
  D(2, 2) = 0.5 - nu;
  D *= E / (1.0 + nu) / (1.0 - 2.0 * nu);

  // Constructing the stiffness matrix
  Mat<double> q4(4, 2);
  q4(0, 0) = q4(0, 1) = q4(1, 1) = q4(2, 0) = -1.0;
  q4(1, 0) = q4(2, 1) = q4(3, 0) = q4(3, 1) = 1.0;
  q4 *= 1.0 / sqrt(3.0);

  cout << "Assembling stifness matrix..." << endl;

  // initilizing the B matrix
  Mat<double> B(3, 8, fill::zeros);

  // loop over elements constructing global stiffness matrix K
  for (int i = 0; i < mesh.num_elements(); i++) {
    Mat<double> xIe(4, 2, fill::zeros);
    Mat<double> curr_elem = mesh.conns_set.row(i); // problematic with vector

    // loop over current element nodes number (Implementation: C++11)
    // storing coordinates in a 4x2 matrix: xIe
    int j = 0;
    for (const auto &val : curr_elem) {
      xIe.row(j) = mesh.nodes_set.row(val);
      j++;
    }

    // initializing the stiffness matrix of current element
    Mat<double> Ke(8, 8, fill::zeros);
    Mat<double> dN(2, 4, fill::zeros);
    for (int i = 0; i < int(q4.n_rows); i++) {
      vec q = trans(q4.row(i)); // has to be transposed!!
      dN = gradshape(q);
      Mat<double> J = dN * xIe;
      dN = inv(J) * dN;
      B = Set_B(dN);

      Ke += det(J) * ((B.t() * D) * B); // current element stiffness matrix
    }

    // Now scatter the Ke into K!
    int k = 0;
    for (const auto &val1 : curr_elem) {
      int j = 0;
      for (const auto &val2 : curr_elem) {
        K(2 * val1, 2 * val2) += Ke(2 * k, 2 * j);
        K(2 * val1 + 1, 2 * val2) += Ke(2 * k + 1, 2 * j);
        K(2 * val1 + 1, 2 * val2 + 1) += Ke(2 * k + 1, 2 * j + 1);
        K(2 * val1, 2 * val2 + 1) += Ke(2 * k, 2 * j + 1);
        j++;
      }
      k++;
    }
  }
}

/*
** Shape function for a 4-node, isoparametric element!
*/
vec SolidElement::shape(vec xi) {
  double x = xi(0);
  double y = xi(1);

  // C++11 implementation
  vec N = {{(1.0 - x) * (1.0 - y), (1.0 + x) * (1.0 - y), (1.0 + x) * (1.0 + y),
            (1.0 - x) * (1.0 + y)}};
  N *= 0.25;
  return N;
}

/*
** Gradient of the shape function for a 4-node, isoparametric element!
*/
mat SolidElement::gradshape(vec xi) {
  double x = xi(0);
  double y = xi(1);

  // C++11 implementation
  Mat<double> dN = {{-(1.0 - y), (1.0 - y), (1.0 + y), -(1.0 + y)},
                    {-(1.0 - x), -(1.0 + x), (1.0 + x), (1.0 - x)}};
  dN *= 0.25;
  return dN;
}

/*
** Set B matrix for a 2D element
*/
mat SolidElement::Set_B(mat dN) {
  Mat<double> B(3, 8, fill::zeros);
  B(0, 0) = dN(0, 0);
  B(0, 2) = dN(0, 1);
  B(0, 4) = dN(0, 2);
  B(0, 6) = dN(0, 3);

  B(1, 1) = dN(1, 0);
  B(1, 3) = dN(1, 1);
  B(1, 5) = dN(1, 2);
  B(1, 7) = dN(1, 3);

  B(2, 0) = dN(1, 0);
  B(2, 1) = dN(0, 0);
  B(2, 2) = dN(1, 1);
  B(2, 3) = dN(0, 1);
  B(2, 4) = dN(1, 2);
  B(2, 5) = dN(0, 2);
  B(2, 6) = dN(1, 3);
  B(2, 7) = dN(0, 3);

  return B;
}
