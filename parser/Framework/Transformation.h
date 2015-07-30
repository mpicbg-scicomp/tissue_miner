#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

#include <stdlib.h>
#include <math.h>
#include "Object2D.h"
#include "Matrix2x2.h"

// class AbstractTransformation {
// public:
//   AbstractTransformation() : _inverse(NULL) {}
//   virtual ~AbstractTransformation() { if(_inverse) delete _inverse; }
//   
//   virtual double mapX(const double x, const double y) const { return map(Vector2D(x, y)).x(); }
//   virtual double mapY(const double x, const double y) const { return map(Vector2D(x, y)).y(); }
//   virtual double invertX(const double x, const double y) const { return invert(Vector2D(x, y)).x(); }
//   virtual double invertY(const double x, const double y) const { return invert(Vector2D(x, y)).y(); }
//   
//   virtual Vector2D map(const Vector2D &v) const = 0;
//   virtual Vector2D invert(const Vector2D &v) const = 0;
//   const AbstractTransformation *getInverse() const {
//     if(!_inverse) { _inverse = generateInverse(); }
//     return _inverse;
//   }
//   
// protected:  
//   virtual const AbstractTransformation *generateInverse() const = 0;
//   
// private:
//   mutable const AbstractTransformation *_inverse;
// };


// map: v -> mv+v0
class MatrixTransformation {
public:
  MatrixTransformation(const Vector2D &r0, const Matrix2x2 &m) : _x0(r0.x()), _y0(r0.y()), _m11(m.m11()), _m12(m.m12()), _m21(m.m21()), _m22(m.m22()), _inverse(NULL) {}
  MatrixTransformation(double x0, double y0, double m11, double m12, double m21, double m22) : _x0(x0), _y0(y0), _m11(m11), _m12(m12), _m21(m21), _m22(m22), _inverse(NULL) {}
  MatrixTransformation() : _x0(0), _y0(0), _m11(1), _m12(0), _m21(0), _m22(1), _inverse(NULL) {}
  MatrixTransformation(const MatrixTransformation &other) : _x0(other._x0), _y0(other._y0), _m11(other._m11), _m12(other._m12), _m21(other._m21), _m22(other._m22), _inverse(NULL) {}
  ~MatrixTransformation() { if(_inverse) { delete _inverse; } }
  
  double mapX(const double x, const double y) const { return _m11*x+_m12*y+_x0; }
  double mapY(const double x, const double y) const { return _m21*x+_m22*y+_y0; }
  Vector2D map(const Vector2D &v) const { return Vector2D(mapX(v.x(), v.y()), mapY(v.x(), v.y())); }
  Vector2D mapLocally(const Vector2D &v) const { return Vector2D(_m11*v.x()+_m12*v.y(), _m21*v.x()+_m22*v.y()); }
  void map(Vector2D &vT, const Vector2D &v) const { vT.set(mapX(v.x(), v.y()), mapY(v.x(), v.y())); }
  void mapLocally(Vector2D &vT, const Vector2D &v) const { vT.set(_m11*v.x()+_m12*v.y(), _m21*v.x()+_m22*v.y()); }

  double invertX(const double x, const double y) const { return getInverse()->mapX(x, y); }
  double invertY(const double x, const double y) const { return getInverse()->mapY(x, y); }
  Vector2D invert(const Vector2D &v) const { return getInverse()->map(v); }
  Vector2D invertLocally(const Vector2D &v) const { return getInverse()->mapLocally(v); }
  void invert(Vector2D &vT, const Vector2D &v) const { getInverse()->map(vT, v); }
  void invertLocally(Vector2D &vT, const Vector2D &v) const { getInverse()->mapLocally(vT, v); }

  double dxt_dx(const double x, const double y) const { return _m11; }
  double dxt_dy(const double x, const double y) const { return _m12; }
  double dyt_dx(const double x, const double y) const { return _m21; }
  double dyt_dy(const double x, const double y) const { return _m22; }

  // (a*b).map(v) == a.map(b.map(v))
  friend MatrixTransformation operator*(const MatrixTransformation &lhs, const MatrixTransformation &rhs) { 
    return MatrixTransformation(lhs.mapX(rhs._x0, rhs._y0), lhs.mapY(rhs._x0, rhs._y0), lhs._m11*rhs._m11+lhs._m12*rhs._m21, lhs._m11*rhs._m12+lhs._m12*rhs._m22, lhs._m21*rhs._m11+lhs._m22*rhs._m21, lhs._m21*rhs._m12+lhs._m22*rhs._m22);
  }

  virtual const MatrixTransformation *inverse() const { return _inverse; }
  const MatrixTransformation *getInverse() const {
    if(!_inverse) { _inverse = generateInverse(); }
    return _inverse;
  }

protected:  
  virtual const MatrixTransformation *generateInverse() const {
    double recOfDet = 1.0/(_m11*_m22 - _m12*_m21);
    return new MatrixTransformation((_m12*_y0-_m22*_x0)*recOfDet, (_m21*_x0-_m11*_y0)*recOfDet, _m22*recOfDet, -_m12*recOfDet, -_m21*recOfDet, _m11*recOfDet);
  }
  
private:
  double _x0, _y0;
  double _m11, _m12, _m21, _m22;
  mutable const MatrixTransformation *_inverse;
};


class IsotropicTransformation : public MatrixTransformation {
public:
  IsotropicTransformation(const Vector2D &deltaR=Vector2D(0.0, 0.0), const double scale=1.0, const double angleInRad=0.0, const bool flip=false) : MatrixTransformation(deltaR.x(), deltaR.y(), scale*cos(angleInRad), -scale*sin(angleInRad)*(flip?-1:1), scale*sin(angleInRad), scale*cos(angleInRad)*(flip?-1:1)), _deltaR(deltaR), _scale(scale), _angle(angleInRad), _flip(flip?-1:1) {}
  
  // (a*b).map(v) == a.map(b.map(v))
  friend IsotropicTransformation operator*(const IsotropicTransformation &lhs, const IsotropicTransformation &rhs) { 
    return IsotropicTransformation(lhs.map(rhs._deltaR), lhs._scale*rhs._scale, lhs._angle+lhs._flip*rhs._angle, lhs._flip*rhs._flip<0);
  }

  virtual IsotropicTransformation *inverse() const { return (IsotropicTransformation*)getInverse();  }

  double angle() const { return _angle; }
  double scale() const { return _scale; }
  double flip() const { return _flip; }
  
  double mapAngle(const double phi) const { return _angle + _flip*phi; }
  double invertAngle(const double phiT) const { return inverse()->mapAngle(phiT); }

protected:
  virtual const IsotropicTransformation *generateInverse() const {
    return new IsotropicTransformation(IsotropicTransformation(Vector2D(0,0), -1.0/_scale, -_flip*_angle, _flip<0).map(_deltaR), 1.0/_scale, -_flip*_angle, _flip<0);
  }
  
private:
  Vector2D _deltaR;
  double _scale, _angle, _flip;
};


class Rotation : public IsotropicTransformation {
public:
  Rotation(const Vector2D originInOldCoords, const double angleInRad=0.0, const bool flip=false) : 
    IsotropicTransformation(-rotatedObject(originInOldCoords, angleInRad), 1.0, angleInRad, flip) {}
};


// class IdentityTransformation : public Rotation {
// public:
//   virtual Vector2D map(const Vector2D &v) { return v; }
//   virtual Vector2D invert(const Vector2D &v) { return v; }
//   virtual double dxt_dx(const double x, const double y) const { return 1; }
//   virtual double dxt_dy(const double x, const double y) const { return 0; }
//   virtual double dyt_dx(const double x, const double y) const { return 0; }
//   virtual double dyt_dy(const double x, const double y) const { return 1; }
// protected:
//   virtual const IdentityTransformation *generateInverse() const { return new IdentityTransformation; }
// };

extern IsotropicTransformation Identity;

#endif
