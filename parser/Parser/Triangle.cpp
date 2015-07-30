#include <math.h>
#include "Matrix2x2.h"
#include "Triangle.h"

const double ReferenceSide0 = 2.0/sqrt(sqrt(3)); // side of a regular triangle with area 1
const Matrix2x2 ReferenceTriangleMatrix = ReferenceSide0*Matrix2x2(-0.5*sqrt(3), -0.5*sqrt(3), 0.5, -0.5); // the columns of this matrix are the two vectors r1-r0 and r2-r0; where ri are the corners of a regular triangle, sorted in ccw direction; the area of the triangle is 1.
const Matrix2x2 ReferenceTriangleInverseMatrix = ReferenceTriangleMatrix.inverse();

Triangle::Triangle(const Vector2D &r0, const Vector2D &r1, const Vector2D &r2, const Triangle *const dT) : t(Matrix2x2(r1.x()-r0.x(), r2.x()-r0.x(), r1.y()-r0.y(), r2.y()-r0.y())*ReferenceTriangleInverseMatrix), area(t.determinant()), deformedTriangle(dT) {
  if(!t.decomposeIntoScalingRotationShear(s, theta, Q)) {
//       cout << "Triangle::Triangle: Zero determinant!" << endl;
  }
}
