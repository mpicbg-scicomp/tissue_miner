#ifndef GRID_IN_SPACE_H
#define GRID_IN_SPACE_H

#include "Object2D.h"
#include "Transformation.h"

class GridInSpace {
public:
  GridInSpace(double minX, double maxX, double minY, double maxY, double boxSizeX, double boxSizeY) : _minX(minX), _minY(minY), _boxSizeX(boxSizeX), _boxSizeY(boxSizeY) {
    _numBoxesX = ceil((maxX-_minX)/_boxSizeX);
    _numBoxesY = ceil((maxY-_minY)/_boxSizeY);
  }

  // integer grid:
  int numBoxesX() const { return _numBoxesX; }
  int numBoxesY() const { return _numBoxesY; }
  bool checkSize(const int i, const int j) const {
    return (i>=0) && (i<numBoxesX()) && (j>=0) && (j<numBoxesY());
  }

  // grid in space:
  double boxSizeX() const { return _boxSizeX; }
  double boxSizeY() const { return _boxSizeY; }
  double minBoxSize() const { return (_boxSizeX<_boxSizeY)?_boxSizeX:_boxSizeY; }
  
  // conversion of positions
  Vector2D boxMid(const int column, const int row) const { return Vector2D((column+0.5)*_boxSizeX+_minX, (row+0.5)*_boxSizeY+_minY); }
  Vector2D boxCornerBottomLeft(const int column, const int row) const { return Vector2D((column)*_boxSizeX+_minX, (row)*_boxSizeY+_minY); }
  Vector2D boxCornerBottomRight(const int column, const int row) const { return Vector2D((column+1)*_boxSizeX+_minX, (row)*_boxSizeY+_minY); }
  Vector2D boxCornerTopLeft(const int column, const int row) const { return Vector2D((column)*_boxSizeX+_minX, (row+1)*_boxSizeY+_minY); }
  Vector2D boxCornerTopRight(const int column, const int row) const { return Vector2D((column+1)*_boxSizeX+_minX, (row+1)*_boxSizeY+_minY); }
  void boxIndex(int &column, int &row, const Vector2D v) const { 
    column = (v.x()-_minX)/_boxSizeX;
    row = (v.y()-_minY)/_boxSizeY;
  }
  bool checkSize(const Vector2D v) const {
    int column, row;
    boxIndex(column, row, v);
    return checkSize(column, row);
  }
  
  // coarsen
  GridInSpace coarsen(int nx, int ny=-1) const {
    if(ny<0) { ny = nx; }
    return GridInSpace(_minX, _minX+numBoxesX()*boxSizeX(), _minY, _minY+numBoxesY()*boxSizeY(), nx*boxSizeX(), ny*boxSizeY());
  }
  
private:
  double _minX, _minY;
  int _numBoxesX, _numBoxesY;
  double _boxSizeX, _boxSizeY;
};

#endif
