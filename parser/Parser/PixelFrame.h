/* 
 * File:   PixelFrame.h
 * Author: mmpi
 *
 * Created on September 13, 2013, 6:42 PM
 */

#ifndef PIXELFRAME_H
#define	PIXELFRAME_H

#include <sstream>
#include <string>
#include "Object2D.h"
#include "AbstractRasterCanvas.h"

/** This class provides basic routines to quickly handle pixels inside of a raster image */
class PixelFrame {
public:
  /**
   * creates a pixel frame from a raster canvas
   * @param canvas raster canvas, may be destroyed afterwards
   */
  PixelFrame(const AbstractRasterCanvas &canvas) : Width(canvas.width()), Height(canvas.height()), TotalPixelNumber(canvas.width()*canvas.height()) {
    for(int i = 0; i < NumberOfNeighbors; ++i) {
      _neighborOffset[i] = NeighborOffsetX[i]*Height + NeighborOffsetY[i];
    }
    _imageData = canvas.rgbData();
  }
  ~PixelFrame() { delete[] _imageData; }

  /** is pixel in image anymore? */
  bool inCanvas(const int p) const { return p<TotalPixelNumber; }
  /** fraction of the whole canvas? */
  double fractionOfCanvas(const int p) const { return 1.0*p/TotalPixelNumber; }
  /** is pixel on the image margin (i.e. on the outermost pixel rows or columns)? */
  bool isOnMargin(const int p) const { 
    return (p<Height) || (p>=TotalPixelNumber-Height) || (p%Height==0) || (p%Height==Height-1);
  }
  /** returns the neighbor offsets in ccw direction, starting with "right" */
  int isNeighborOnCanvas(const int p, const int index) const {
    return ((p>=Height) || (NeighborOffsetX[index]>=0))
            && ((p<TotalPixelNumber-Height) || (NeighborOffsetX[index]<=0))
            && ((p%Height>0) || (NeighborOffsetY[index]>=0))
            && ((p%Height<Height-1) || (NeighborOffsetY[index]<=0));
  }

  /**
   * @param p pixel
   * @return pixel position string
   */
  std::string positionString(const int p) const {
    std::stringstream ss;
//    ss << "(" << p/Height << ", " << p%Height << ")";
    ss << "(" << p/Height << ", " << Height-1-p%Height << ")"; // for debugging
    return ss.str();
  }
  /**
   * @param p pixel
   * @return 2D vector corresponding to pixel p
   */
  Vector2D toVector(const int p) const {
    return Vector2D(p/Height, p%Height);
  }
  int toPosition(const Vector2D &v) const {
    return int(v.x()+0.5)*Height + int(v.y()+0.5);
  }
  
  /** the number of neighbor pixels */
  static const int NumberOfNeighbors = 8;
  /** returns the neighbor offsets in ccw direction, starting with "right" */
  int direction(const int index) const { return _neighborOffset[index]; }
  /** returns true, if the direction is diagonal */
  static bool diagonalNeighbor(const int index) { return index & 1; }
  /** returns the opposite direction */
  static int oppositeDirection(const int index) { return (index + 4) % NumberOfNeighbors; }
  /** returns distance to neighbor pixel */
  static const double NNeighborDistance[NumberOfNeighbors];
  
  /** gets image data (pixel color in format 0x00RRGGBB) */
  unsigned int data(const int p) const { return _imageData[p]; }
  
  
private:
  const int Width, Height, TotalPixelNumber;
  static const int NeighborOffsetX[NumberOfNeighbors];
  static const int NeighborOffsetY[NumberOfNeighbors];
  int _neighborOffset[NumberOfNeighbors];
  unsigned int *_imageData;
};

#endif	/* PIXELFRAME_H */

