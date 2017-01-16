#ifndef _SOLIDELEMENT_H_
#define _SOLIDELEMENT_H_

#include "Mesh.h"
#include <armadillo>

class SolidElement {
public:
  SolidElement(Mesh, double, double);
  virtual ~SolidElement();
  mat Set_B(mat);
  vec shape(vec);
  mat gradshape(vec);

  mat K;
};

#endif
