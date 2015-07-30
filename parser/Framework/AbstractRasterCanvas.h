#ifndef ABSTRACT_RASTER_CANVAS_H
#define ABSTRACT_RASTER_CANVAS_H

#include <stdlib.h>

#include "Color.h"
#include "AbstractCanvas.h"
// #include "ComplexImage.h"

// j=0 is at the bottom...
class AbstractRasterCanvas : public virtual AbstractCanvas {
public:
  virtual bool load(const char filename[]) = 0;
  virtual bool save(const char filename[]) const = 0;

  virtual void drawLineD(const double x1, const double y1, const double x2, const double y2, const Color &col);
  virtual void drawPixel(const int i, const int j, const Color &col) = 0;
  
  // line-wise; caller has to care for delete[] !!!
  virtual unsigned char* redChannel() const = 0;
  virtual unsigned char* greenChannel() const = 0;
  virtual unsigned char* blueChannel() const = 0;
  virtual unsigned char* alphaChannel() const { return NULL; };
  virtual unsigned int* rgbaData() const = 0;
  virtual unsigned int* rgbData() const = 0;
  unsigned char* channel(const ColorChannel c) const;
//   ComplexImage channelToComplexImage(const ColorChannel c) const {
//     unsigned char *data = channel(c);
//     ComplexImage cI(width(), height(), data);
//     delete[] data;
//     return cI;
//   };
};

#endif
