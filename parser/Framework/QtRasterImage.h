#ifndef QT_RASTER_IMAGE_H
#define QT_RASTER_IMAGE_H

#include <Qt/QtCore>
#include <Qt/QtGui>

#include "Color.h"
#include "Transformation.h"
#include "QtImage.h"
#include "AbstractRasterCanvas.h"

// j=0 is at the bottom...
class QtRasterImage : public QtImage, public AbstractRasterCanvas {
public:
  QtRasterImage();
  QtRasterImage(const char filename[]);
  QtRasterImage(const int width, const int height);
  QtRasterImage(const QtRasterImage &img);
  ~QtRasterImage();

  bool valid() const { return _valid; }
  virtual void reset(const Color &col=Black);
  void invert();

  virtual bool load(const char filename[]);
  virtual bool save(const char filename[]) const;

  void setRedChannel(const unsigned char r);
  void setGreenChannel(const unsigned char g);
  void setBlueChannel(const unsigned char b);
  void setAlphaChannel(const unsigned char a);

  void copyToRedChannel(const unsigned char *const data);
  void copyToGreenChannel(const unsigned char *const data);
  void copyToBlueChannel(const unsigned char *const data);
  void copyToAlphaChannel(const unsigned char *const data);

  // line-wise; caller has to care for delete[] !!!
  virtual unsigned char* redChannel() const;
  virtual unsigned char* greenChannel() const;
  virtual unsigned char* blueChannel() const;
  virtual unsigned char* alphaChannel() const;
  unsigned int* rgbaData() const;
  unsigned int* rgbData() const;

  bool inImage(const int i, const int j) const { return (i>=0) && (i<width()) && (j>=0) && (j<height()); }
  char pixel(const int i, const int j, const ColorChannel c) const;
  Color pixel(const int i, const int j) const {
    return Color(pixel(i, j, RedChannel), pixel(i, j, GreenChannel), pixel(i, j, BlueChannel), pixel(i, j, AlphaChannel));
  }
  virtual void drawLineD(const double x1, const double y1, const double x2, const double y2, const Color &col) {
    QtImage::drawLineD(x1, y1, x2, y2, col);
  }
  virtual void drawPixel(const int i, const int j, const Color &col);
  void blend(const MatrixTransformation &t, const QtRasterImage &img, const ColorChannel targetChannel=NoChannel, const ColorChannel sourceChannel=NoChannel);
  QtRasterImage scaled(int width, int height);
  
protected:
  virtual double yTransformation(const double y) const;

  bool _valid;
  QImage *_image;
  QPainter *_painter;
  virtual QPaintDevice *paintDevice() const { return _image; }
  virtual QPainter *startPainting() { return _painter = new QPainter(_image); }
  virtual void endPainting() { delete _painter; }

  friend class QtVectorImage;
};

#endif
