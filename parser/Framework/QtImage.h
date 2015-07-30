#ifndef QT_IMAGE_H
#define QT_IMAGE_H

#include <vector>
#include <Qt/QtCore>
#include <Qt/QtGui>

#include "Color.h"
#include "AbstractCanvas.h"
#include "AbstractField.h"

class MyQColor : public QColor {
public:
  MyQColor(const Color &col) : QColor(col.r(), col.g(), col.b(), col.a()) {};
};

class QtImage : public virtual AbstractCanvas {
public:
  QtImage() : _antialiasing(true) {}
  virtual int width() const { return paintDevice()->width(); }
  virtual int height() const { return paintDevice()->height(); }
  void setAntialiasing(bool on=true) { _antialiasing=on; }

  virtual void drawLineD(const double x1, const double y1, const double x2, const double y2, const Color &col);
  virtual void drawBarD(const double x1, const double y1, const double x2, const double y2, const Color &col, const double size=2.0);
  virtual void drawArrowD(const double midX, const double midY, const double deltaX, const double deltaY, const Color &col, const double size=2.5);

  void drawPolygon(const std::vector<Vector2D> &vertices, const Color &col);
  void drawPolygonWithOutline(const std::vector<Vector2D> &vertices, const Color &col, const double width, const Color &outerCol);
  void drawRectD(const double x1, const double y1, const double w, const double h, const Color &col);
  void drawArrowExt(const double midX, const double midY, const double deltaX, const double deltaY, const Color &col, const Color &headCol, const double size=2.5);
  void drawArrowExtC(const double midX, const double midY, const double deltaX, const double deltaY, const Color &col, const Color &headCol, const Color &cCol, const double size=2.5);
  void drawText(const double x, const double y, const double pixelSize, const Color &col, const char* text);
  void drawCircle(const double x, const double y, const double radius, const Color &col);
  void drawEmptyCircle(const double x, const double y, const double radius, const Color &col, const double thickness=1.0);
  void drawEllipse(const double x, const double y, const double radius1, const double radius2, const double angle, const Color &col);

  void drawRect(const Vector2D &lowerLeft, const double w, const double h, const Color &col) { drawRectD(lowerLeft.x(), lowerLeft.y(), w, h, col); }
  void drawArrowExt(const Vector2D &mid, const Vector2D &delta, const Color &col, const Color &headCol, const double size=2.5) { drawArrowExt(mid.x(), mid.y(), delta.x(), delta.y(), col, headCol, size); }
  void drawArrowExtC(const Vector2D &mid, const Vector2D &delta, const Color &col, const Color &headCol, const Color &cCol, const double size=2.5) { drawArrowExtC(mid.x(), mid.y(), delta.x(), delta.y(), col, headCol, cCol, size); }
  void drawText(const Vector2D &lowerLeft, const double pixelSize, const Color &col, const char* text) { drawText(lowerLeft.x(), lowerLeft.y(), pixelSize, col, text); }
  void drawCircle(const Vector2D &mid, const double radius, const Color &col) { drawCircle(mid.x(), mid.y(), radius, col); }
  void drawEmptyCircle(const Vector2D &mid, const double radius, const Color &col, const double thickness=1.0) { drawEmptyCircle(mid.x(), mid.y(), radius, col, thickness); }
  void drawEllipse(const Vector2D &mid, const double radius1, const double radius2, const double angle, const Color &col) { drawEllipse(mid.x(), mid.y(), radius1, radius2, angle, col); }

  void drawScalarArray(const MatrixTransformation &trafo, const AbstractScalarArray &a, const ColorBar &Cb=Rainbow, const unsigned char AlphaValue=0x80);
  void drawVector2DArray(const IsotropicTransformation &trafo, const AbstractVector2DArray &a, const bool normalized=true, const Color &Col=Blue, const Color &HeadCol=Red, const double RelLength=0.9, const double RelWidth=0.35);
  
protected:
  virtual double yTransformation(const double y) const = 0;

  bool _antialiasing;
  virtual QPaintDevice *paintDevice() const = 0;
  virtual QPainter *startPainting() = 0;
  virtual void endPainting() = 0;
};

#endif
