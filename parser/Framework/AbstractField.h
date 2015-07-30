#ifndef ABSTRACT_FIELD_H
#define ABSTRACT_FIELD_H

#include <iostream>
#include "Object2D.h"
#include "GridInSpace.h"

template<class T> class AbstractField {
public:
  virtual bool valueAt(T &v, const Vector2D &r) const = 0;
  virtual bool valid(const Vector2D &r) const = 0;
};

typedef AbstractField<Vector2D> AbstractVector2DField;
typedef AbstractField<Nematic2D> AbstractNematic2DField;

template<class T> class AbstractArray : public AbstractField<T> {
public:
  AbstractArray(const GridInSpace &grid) : _grid(grid) { }
  const GridInSpace &grid() const { return _grid; }

  // may be overloaded:
  virtual bool valid(const int i, const int j) const { return _grid.checkSize(i, j); }
  virtual T operator()(const int i, const int j) const = 0;
  
  virtual bool valueAt(T &v, const Vector2D &r) const {
    int i, j;
    _grid.boxIndex(i, j, r);
    if(valid(i,j)) {
      v = (*this)(i,j);
      return true;
    } else {
      return false;
    }
  }
  virtual bool valid(const Vector2D &r) const { return _grid.checkSize(r); };
  
private:
  const GridInSpace &_grid;
};

typedef AbstractArray<double> AbstractScalarArray;
typedef AbstractArray<Vector2D> AbstractVector2DArray;
typedef AbstractArray<Nematic2D> AbstractNematic2DArray;

// constant value

template<class T> class ConstantArray : public AbstractArray<T> {
public:
  ConstantArray(const T value, const GridInSpace &grid) : AbstractArray<T>(grid), _value(value) { }
  virtual T operator()(int i, int j) const { return _value; }
private:
  const T _value;
};

typedef ConstantArray<double> ConstantScalarArray;
typedef ConstantArray<Vector2D> ConstantVector2DArray;
typedef ConstantArray<Nematic2D> ConstantNematic2DArray;


// operations

template<class T> class NormalizedArray : public AbstractArray<T> {
public:
  NormalizedArray(const AbstractArray<T> &arg) : AbstractArray<T>(arg.grid()), _arg(arg) { }
  virtual T operator()(int i, int j) const { return _arg(i,j)/_arg(i,j).norm(); }
  virtual bool valid(int i, int j) const { return _arg.valid(i,j); }
private:
  const AbstractArray<T> &_arg;
};

typedef NormalizedArray<Vector2D> NormalizedVector2DArray;
typedef NormalizedArray<Nematic2D> NormalizedNematic2DArray;

template<class T> class NormalizedSaturatedArray : public AbstractArray<T> {
public:
  NormalizedSaturatedArray(const AbstractArray<T> &other, const double normalization) : AbstractArray<T>(other.grid()), _other(other), _normalization(normalization) {}
  virtual T operator()(int i, int j) const {
    T b(_other(i,j));
    if(b.norm()>_normalization) {
      return b/b.norm();
    } else {
      return b/_normalization;
    }
  }
  virtual bool valid(int i, int j) const { return _other.valid(i,j); }
private:
  const AbstractArray<T> &_other;
  const double _normalization;
};

template<class T, class U, class V> class ProductArray : public AbstractArray<T> {
public:
  ProductArray(const AbstractArray<U> &lhs, const AbstractArray<V> &rhs) : AbstractArray<T>(lhs.grid()), _lhs(lhs), _rhs(rhs) {
    if(&lhs.grid()!=&rhs.grid()) {
      std::cout << "Different grids in ProductArray!" << std::endl;
      exit(1);
    }
  }
  virtual T operator()(int i, int j) const { return _lhs(i,j)*_rhs(i,j); }
  virtual bool valid(int i, int j) const { return _lhs.valid(i,j) && _rhs.valid(i,j); }
private:
  const AbstractArray<U> &_lhs;
  const AbstractArray<V> &_rhs;
};

// template<class T, class U, class V> inline ProductArray<T,U,V> operator*(const AbstractArray<U> &lhs, const AbstractArray<V> &rhs) {
//   return ProductArray<T,U,V>(lhs, rhs);
// }

inline ProductArray<double,double,double> operator*(const AbstractArray<double> &lhs, const AbstractArray<double> &rhs) {
  return ProductArray<double,double,double>(lhs, rhs);
}
inline ProductArray<double,Vector2D,Vector2D> operator*(const AbstractArray<Vector2D> &lhs, const AbstractArray<Vector2D> &rhs) {
  return ProductArray<double,Vector2D,Vector2D>(lhs, rhs);
}
inline ProductArray<double,Nematic2D,Nematic2D> operator*(const AbstractArray<Nematic2D> &lhs, const AbstractArray<Nematic2D> &rhs) {
  return ProductArray<double,Nematic2D,Nematic2D>(lhs, rhs);
}

inline ProductArray<Vector2D,Vector2D,Nematic2D> operator*(const AbstractArray<Vector2D> &lhs, const AbstractArray<Nematic2D> &rhs) {
  return ProductArray<Vector2D,Vector2D,Nematic2D>(lhs, rhs);
}
inline ProductArray<Vector2D,Nematic2D,Vector2D> operator*(const AbstractArray<Nematic2D> &lhs, const AbstractArray<Vector2D> &rhs) {
  return ProductArray<Vector2D,Nematic2D,Vector2D>(lhs, rhs);
}


#endif
