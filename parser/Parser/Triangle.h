/* 
 * File:   Triangle.h
 * Author: mmpi
 *
 * Created on 30. September 2013, 18:54
 */

#ifndef TRIANGLE_H
#define	TRIANGLE_H

#include <stdlib.h>
#include "Matrix2x2.h"

struct Triangle {
  Triangle(const Vector2D &r0, const Vector2D &r1, const Vector2D &r2, const Triangle *const dT=NULL);
  const Matrix2x2 t;
  const double area;
  double s, theta;
  Nematic2D Q;
  const Triangle *const deformedTriangle;
};

#endif	/* TRIANGLE_H */

