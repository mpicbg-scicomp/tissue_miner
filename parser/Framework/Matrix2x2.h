#ifndef MATRIX_2X2_H
#define MATRIX_2X2_H

#include <stdio.h>
#include <math.h>
#include "Object2D.h"

class Matrix2x2 {
public:
  Matrix2x2() : _m11(1.0), _m12(0.0), _m21(0.0), _m22(1.0) {}
//   Matrix2x2(const double trace) : _m11(trace), _m12(0.0), _m21(0.0), _m22(trace) {}
  Matrix2x2(const double m11, const double m12, const double m21, const double m22) : _m11(m11), _m12(m12), _m21(m21), _m22(m22) {}
  Matrix2x2(const Matrix2x2 &other) : _m11(other._m11), _m12(other._m12), _m21(other._m21), _m22(other._m22) {}
  Matrix2x2(const Nematic2D &n) : _m11(n.c1()), _m12(n.c2()), _m21(n.c2()), _m22(-n.c1()) {}
  
  // read
  double m11() const { return _m11; }
  double m12() const { return _m12; }
  double m21() const { return _m21; }
  double m22() const { return _m22; }
  
  // components
  double trace() const { return (_m11+_m22); }
  Matrix2x2 tracePart() const { double tr=0.5*trace(); return Matrix2x2(tr, 0.0, 0.0, tr); }
  Nematic2D symmetricTracelessPart() const { return Nematic2D(0.5*(_m11-_m22), 0.5*(_m12+_m21)); }
  double asymmetry() const { return 0.5*(_m12-_m21); }
  Matrix2x2 asymmetricPart() const { double as=asymmetry(); return Matrix2x2(0.0, as, -as, 0.0); }
  bool decomposeIntoScalingRotationShear(double &s, double &theta, Nematic2D &Q) const;
  
  // inverses
  Matrix2x2 operator-() const { Matrix2x2 tmp(*this); tmp*=-1; return tmp; }
  double determinant() const { return _m11*_m22 - _m12*_m21; }
  Matrix2x2 inverse() const { return Matrix2x2(_m22, -_m12, -_m21, _m11)/determinant(); }

  // manipulate
  void set(const double m11, const double m12, const double m21, const double m22) { _m11 = m11; _m12 = m12; _m21 = m21; _m22 = m22; }
  void operator+=(const Matrix2x2 &o) { _m11+=o._m11; _m12+=o._m12; _m21+=o._m21; _m22+=o._m22; }
  void operator-=(const Matrix2x2 &o) { _m11-=o._m11; _m12-=o._m12; _m21-=o._m21; _m22-=o._m22; }
  void operator*=(const double d) { _m11*=d; _m12*=d; _m21*=d; _m22*=d; }
  void operator*=(const Matrix2x2 &o) { set(_m11*o._m11 + _m12*o._m21, 
                                            _m11*o._m12 + _m12*o._m22, 
                                            _m21*o._m11 + _m22*o._m21, 
                                            _m21*o._m12 + _m22*o._m22); }
  void operator/=(const double d) { _m11/=d; _m12/=d; _m21/=d; _m22/=d; }
  void operator/=(const Matrix2x2 &o) { *this *= o.inverse(); }

  // operations
  Matrix2x2 transposed() const { return Matrix2x2(_m11, _m21, _m12, _m22); }
  friend Matrix2x2 operator+(const Matrix2x2 &lhs, const Matrix2x2 &rhs) { Matrix2x2 tmp(lhs); tmp+=rhs; return tmp; }
  friend Matrix2x2 operator-(const Matrix2x2 &lhs, const Matrix2x2 &rhs) { Matrix2x2 tmp(lhs); tmp-=rhs; return tmp; }
  friend Matrix2x2 operator*(const Matrix2x2 &lhs, const double d) { Matrix2x2 tmp(lhs); tmp*=d; return tmp; }
  friend Matrix2x2 operator*(const double d, const Matrix2x2 &rhs) { Matrix2x2 tmp(rhs); tmp*=d; return tmp; }
  friend Vector2D operator*(const Matrix2x2 &lhs, const Vector2D &rhs) { return Vector2D(lhs._m11*rhs.x()+lhs._m12*rhs.y(), lhs._m21*rhs.x()+lhs._m22*rhs.y()); }
  friend Vector2D operator*(const Vector2D &lhs, const Matrix2x2 &rhs) { return Vector2D(lhs.x()*rhs._m11+lhs.y()*rhs._m21, lhs.x()*rhs._m12+lhs.y()*rhs._m22); }
  friend Matrix2x2 operator*(const Matrix2x2 &lhs, const Matrix2x2 &rhs) { Matrix2x2 tmp(lhs); tmp*=rhs; return tmp; }
  friend Matrix2x2 operator/(const Matrix2x2 &lhs, const double d) { Matrix2x2 tmp(lhs); tmp/=d; return tmp; }
  friend Matrix2x2 operator/(const Matrix2x2 &lhs, const Matrix2x2 &rhs) { Matrix2x2 tmp(lhs); tmp/=rhs; return tmp; }
  
  // print
  void print(const char *format="%7.4f") const {
    char tmp[1024];
    sprintf(tmp, "[%s %s]\n", format, format);
    printf(tmp, _m11, _m12);
    printf(tmp, _m21, _m22);
    printf("\n");
  }
  
private:
  double _m11, _m12, _m21, _m22;
};

extern const Matrix2x2 Identity2x2;

inline Matrix2x2 dyadicProduct(const Vector2D &lhs, const Vector2D &rhs) { return Matrix2x2(lhs.x()*rhs.x(), lhs.x()*rhs.y(), lhs.y()*rhs.x(), lhs.y()*rhs.y()); }
inline Matrix2x2 rotationMatrix(const double angle) { double c=cos(angle), s=sin(angle); return Matrix2x2(c, -s, s, c); }
inline Matrix2x2 exp(const Nematic2D &n) {
  const double nNorm = n.norm();
  if(nNorm==0) {
    return cosh(nNorm)*Identity2x2;
  } else {
    return cosh(nNorm)*Identity2x2 + sinh(nNorm)/nNorm*n;
  }
}

#endif
