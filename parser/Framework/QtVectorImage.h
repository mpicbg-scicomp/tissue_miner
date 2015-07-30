#ifndef QT_VECTOR_IMAGE_H
#define QT_VECTOR_IMAGE_H

#include <Qt/QtCore>
#include <Qt/QtGui>

#include "Transformation.h"
#include "QtImage.h"
#include "QtRasterImage.h"

// j=0 is at the bottom...
class QtVectorImage : public QtImage {
public:
  QtVectorImage(const QtRasterImage &img);
  QtVectorImage(const int width, const int height);
  ~QtVectorImage();

  bool valid() const { return _valid; }

  virtual void reset(const Color &col=White);
  void drawImage(const QtRasterImage &img, const double x=0.0, const double y=0.0);
  void drawImage(const QtRasterImage &img, const Vector2D &lowerLeft) { drawImage(img, lowerLeft.x(), lowerLeft.y()); }
  bool start(const char filename[]);
  void done();

protected:
  virtual double yTransformation(const double y) const;

  void init(const int width, const int height);
  bool _valid;
  QPrinter *_printer;
  QPainter *_painter;

  virtual QPaintDevice *paintDevice() const { return _printer; }
  virtual QPainter *startPainting() { return _painter; }
  virtual void endPainting() {}
};

#endif
