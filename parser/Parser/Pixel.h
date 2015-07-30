/* 
 * File:   Pixel.h
 * Author: mmpi
 *
 * Created on September 13, 2013, 12:46 PM
 */

#ifndef PIXEL_H
#define	PIXEL_H

#include <set>
#include <vector>
#include "PixelFrame.h"
#include "WrongPixelFrameException.h"

/** cell index type */
typedef unsigned int PixelValue;

/** class for containing and managing pixel positions */
class Pixel {
public:
  /** @param frame should not be deleted */
  Pixel(const PixelFrame &frame) :  _frame(frame), _position(0) {}
  /**
   * @param frame should not be deleted
   * @param position position in image
   */
  Pixel(const PixelFrame &frame, const int position) :  _frame(frame), _position(position) {}

  /** comparison of pixel positions */
  bool operator==(const Pixel &other) const {
    if(&_frame!=&(other._frame)) {
      throw WrongPixelFrameException();
    }
    return _position==other._position;
  }
  /** comparison of pixel positions */
  bool operator!=(const Pixel &other) const {
    if(&_frame!=&(other._frame)) {
      throw WrongPixelFrameException();
    }
    return _position!=other._position;
  }
  /** comparison of pixel positions */
  bool operator<(const Pixel &other) const { 
    if(&_frame!=&(other._frame)) {
      throw WrongPixelFrameException();
    }
    return _position<other._position;
  }

  /** setting of pixel positions */
  Pixel &operator=(const Pixel &other) { 
    if(&_frame!=&(other._frame)) {
      throw WrongPixelFrameException();
    }
    _position=other._position; return *this;
  }
  /** walk through image */
  Pixel operator++() { ++_position; return *this; }
  /** is pixel in image anymore? */
  bool inCanvas() const { return _frame.inCanvas(_position); }
 /** fraction of the whole canvas? */
  double fractionOfCanvas() const { return _frame.fractionOfCanvas(_position); }
  /** is pixel on the image margin (i.e. on the outermost pixel rows or columns)? */
  bool isOnMargin() const { return _frame.isOnMargin(_position); }
  /** returns position string */
  std::string positionString() const { return _frame.positionString(_position); }
  /** 2D vector corresponding to pixel */
  Vector2D toVector() const { return _frame.toVector(_position); }
  
  /** this returns the neighbor pixels in ccw direction, starting with "right" */
  Pixel neighbor(const int index) const { return Pixel(_frame, _position+_frame.direction(index)); }
  /** returns true if neighbor is on canvas */
  bool isNeighborOnCanvas(const int index) const { return _frame.isNeighborOnCanvas(_position, index); }
  /** this returns the neighbor pixels in ccw direction, starting with "right" */
  void goTo(const int index) { _position +=_frame.direction(index); }

  /** gets image data (pixel color in format 0x00RRGGBB) */
  PixelValue data() const { return _frame.data(_position); }
  /** bond color index */
  static const PixelValue BondValue;
  static const PixelValue OutsideCanvasValue;
  static const PixelValue DividingCellValue;
  
  /**
   * Go to the neighbor pixel which is also neighbor of the cw-most pixel with value v. *this becomes the cw-most pixel with that property.
   * @param v
   * @return true if everything is consistent.
   */
//  bool goToLastCommonNeighboringBondPixelInCwOrientation(const PixelValue v);
  /**
   * Returns the neighbor index of the pixel which is also neighbor of the cw-most pixel with value v. Also, returns the set of neighboring cell ids.
   * @param neighbors set of neighboring cell ids
   * @param v
   * @return neighbor index, -1 if nothing found
   */
  int getLastCommonNeighboringBondPixelInCwOrientation(std::vector<PixelValue> &neighbors, const PixelValue v, const int oldDirection) const;
  
  /** recursively builds up set of vertices which represent a thick vertex; implementation assumes that direction 0 is not diagonal */
  void addToThickVertex(std::set<Pixel> &thickVertex);
  /** for debugging, only */
  int position() const { return _position; }
  
private:
  const PixelFrame &_frame;
  int _position;
};

#endif	/* PIXEL_H */

