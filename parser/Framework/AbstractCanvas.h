#ifndef ABSTRACT_CANVAS_H
#define ABSTRACT_CANVAS_H

#include "Color.h"
#include "Object2D.h"
#include "AbstractField.h"

// small y at the bottom...
class AbstractCanvas {
public:
  virtual void reset(const Color &col=Black) = 0;

  virtual int width() const = 0;
  virtual int height() const = 0;

  virtual void drawArrowD(const double midX, const double midY, const double deltaX, const double deltaY, const Color &col, const double size=5.0);
  virtual void drawBarD(const double x1, const double y1, const double x2, const double y2, const Color &col, const double size=2.0);
  virtual void drawLineD(const double x1, const double y1, const double x2, const double y2, const Color &col) = 0;

  virtual void drawArrow(const Vector2D &mid, const Vector2D &delta, const Color &col, const double size=5.0) { drawArrowD(mid.x(), mid.y(), delta.x(), delta.y(), col, size); }
  virtual void drawBar(const Vector2D &r1, const Vector2D &r2, const Color &col, const double size=2.0) { drawBarD(r1.x(), r1.y(), r2.x(), r2.y(), col, size); }
  virtual void drawBar(const Vector2D &mid, const Nematic2D &delta, const Color &col, const double size=2.0) {
    Vector2D v = 0.5*vectorFromAngleAndNorm(delta.angle(), delta.norm());
    drawBar(mid-v, mid+v, col, size);
  }
  virtual void drawLine(const Vector2D &r1, const Vector2D &r2, const Color &col) { drawLineD(r1.x(), r1.y(), r2.x(), r2.y(), col);  }
 
  void drawNematic2DArray(const IsotropicTransformation &trafo, const AbstractNematic2DArray &a, const bool normalized=true, const Color &Col=Green, const double RelativeLength=0.9, const double RelWidth=0.12);
 
  bool isInsideImage(Vector2D v) {
    return (v.x()>=0) && (v.y()>=0) && (v.x()<width()) && (v.y()<height());
  }
};

#endif
