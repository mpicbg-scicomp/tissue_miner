#ifndef OBJECT_2D_H
#define OBJECT_2D_H

#include <ostream>
#include <math.h>
#include <stdio.h>

template<int N> class Object2D {
public:
  Object2D() : _c1(0.0), _c2(0.0) {}
  Object2D(const double c1, const double c2) : _c1(c1), _c2(c2) {}
  Object2D(const Object2D<N> &other) : _c1(other._c1), _c2(other._c2) {}
  
  // read
  double c1() const { return _c1; }
  double c2() const { return _c2; }
  double x() const { return c1(); }
  double y() const { return c2(); }
  double angle() const { return atan2(_c2, _c1)/N; } // in rad
  double normSq() const { return _c1*_c1 + _c2*_c2; }
  double norm() const { return sqrt(normSq()); }
  int order() const { return N; }

  // manipulate
  void set(const double c1, const double c2) { _c1 = c1; _c2 = c2; }
  void rotateBy(const double angleInRad) {
    double c1=_c1, c=cos(N*angleInRad), s=sin(N*angleInRad);
    _c1 = c*c1-s*_c2;
    _c2 = s*c1+c*_c2;
  }
  Object2D<N> operator-() { Object2D<N> tmp(*this); tmp*=-1; return tmp; }
  void operator+=(const Object2D<N> &o) { _c1 += o._c1; _c2 += o._c2; }
  void operator-=(const Object2D<N> &o) { _c1 -= o._c1; _c2 -= o._c2; }
  void operator*=(const double d) { _c1 *= d; _c2 *= d; }
  void operator/=(const double d) { _c1 /= d; _c2 /= d; }

  friend Object2D<N> rotatedObject(const Object2D<N> &lhs, const double angleInRad) { Object2D<N> tmp(lhs); tmp.rotateBy(angleInRad); return tmp; }
  friend Object2D<N> operator+(const Object2D<N> &lhs, const Object2D<N> &rhs) { Object2D<N> tmp(lhs); tmp+=rhs; return tmp; }
  friend Object2D<N> operator-(const Object2D<N> &lhs, const Object2D<N> &rhs) { Object2D<N> tmp(lhs); tmp-=rhs; return tmp; }
  friend Object2D<N> operator-(const Object2D<N> &rhs) { Object2D<N> tmp(-1.0*rhs); return tmp; }
  friend Object2D<N> operator*(const Object2D<N> &lhs, const double d) { Object2D<N> tmp(lhs); tmp*=d; return tmp; }
  friend Object2D<N> operator*(const double d, const Object2D<N> &rhs) { Object2D<N> tmp(rhs); tmp*=d; return tmp; }
  friend double operator*(const Object2D<N> &lhs, const Object2D<N> &rhs) { return lhs.x()*rhs.x() + lhs.y()*rhs.y(); }
  friend Object2D<N> elementWiseProduct(const Object2D<N> &lhs, const Object2D<N> &rhs) { Object2D<N> tmp(lhs.x()*rhs.x(), lhs.y()*rhs.y()); return tmp; }
  friend double crossProductZ(const Object2D<N> &lhs, const Object2D<N> &rhs) { return lhs.x()*rhs.y() - lhs.y()*rhs.x(); }
  friend Object2D<N> operator/(const Object2D<N> &lhs, const double d) { Object2D<N> tmp(lhs); tmp/=d; return tmp; }
  friend std::ostream& operator<<(std::ostream& os, const Object2D<N> &rhs) { return os << "(" << rhs._c1 << ", " << rhs._c2 << ")"; }

  // file IO
  bool readFromFile(FILE *f) {
    double dummy1, dummy2;
    if(2==fscanf(f, "%lg %lg", &dummy1, &dummy2)) { _c1 = dummy1; _c2 = dummy2; return true; }
/*    int dummy1, dummy2;
    if(2==fscanf(f, "%d %d", dummy1, dummy2)) { _c1 = dummy1; _c2 = dummy2; return true; }*/
    return false;
  }
  bool writeToFile(FILE *f) const { return 0<=fprintf(f, "%lg %lg", _c1, _c2); }

private:
  double _c1, _c2;
};

typedef Object2D<1> Vector2D;
typedef Object2D<2> Nematic2D;

template<int N> Object2D<N> objectFromAngleAndNorm(const double angle, const double norm=1.0) { Object2D<N> tmp(norm*cos(N*angle), norm*sin(N*angle)); return tmp; }
inline Vector2D vectorFromAngleAndNorm(const double angle, const double norm=1.0) { return objectFromAngleAndNorm<1>(angle, norm); }
inline Nematic2D nematicFromAngleAndNorm(const double angle, const double norm=1.0) { return objectFromAngleAndNorm<2>(angle, norm); }

inline Vector2D operator*(const Nematic2D &lhs, const Vector2D &rhs) { return Vector2D(lhs.c1()*rhs.x()+lhs.c2()*rhs.y(), lhs.c2()*rhs.x()-lhs.c1()*rhs.y()); }
inline Vector2D operator*(const Vector2D &lhs, const Nematic2D &rhs) { return Vector2D(lhs.x()*rhs.c1()+lhs.y()*rhs.c2(), lhs.x()*rhs.c2()-lhs.y()*rhs.c1()); }
  
#endif
