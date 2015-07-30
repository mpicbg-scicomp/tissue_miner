#include "Matrix2x2.h"

const Matrix2x2 Identity2x2(1.0, 0.0, 0.0, 1.0);

bool Matrix2x2::decomposeIntoScalingRotationShear(double &s, double &theta, Nematic2D &Q) const {
  const double det = determinant();
  s = 0.5*log(fabs(det));
  if(!isfinite(s)) {
    s = theta = 0.0;
    Q.set(0,0);
    return false;
  }
  theta = atan2(-asymmetry(), 0.5*trace());
  const Nematic2D n(symmetricTracelessPart());
  const double twoPhi = theta + 2*n.angle();
  const double normQ = asinh(n.norm()/exp(s));
  Q.set(normQ*cos(twoPhi), normQ*sin(twoPhi));
  return true;
}
