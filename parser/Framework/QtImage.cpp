#include <Qt/QtCore>
#include <Qt/QtGui>

#include "QtImage.h"

void QtImage::drawRectD(const double x1, const double y1, const double w, const double h, const Color &col) {
  QPainter *painter = startPainting();
  painter->setRenderHint(QPainter::Antialiasing, _antialiasing);
  painter->setPen(QPen(Qt::NoPen));
  painter->setBrush(QBrush(MyQColor(col)));
  QRectF rect(x1, yTransformation(h+y1), w, h);
  painter->drawPolygon(rect);
  endPainting();
}


void QtImage::drawBarD(const double x1, const double y1, const double x2, const double y2, const Color &col, const double size) {
  QPainter *painter = startPainting();
  painter->setRenderHint(QPainter::Antialiasing, _antialiasing);
  painter->setPen(QPen(Qt::NoPen));
  painter->setBrush(QBrush(MyQColor(col)));

  double deltaX = x2-x1;
  double deltaY = y2-y1;
  double midX = 0.5*(x1+x2);
  double midY = 0.5*(y1+y2);
  double len = sqrt(deltaX*deltaX + deltaY*deltaY);
 
  double angle = atan2(deltaY, deltaX)/M_PI*180.0;

  if(len<size) {
    painter->drawEllipse(QPointF(midX, yTransformation(midY)), 0.5*size, 0.5*size);
    endPainting();
    return;
  }

  QRectF bar(-0.5*len, -0.5*size, len, size);
  painter->drawPolygon(QMatrix().translate(midX, yTransformation(midY)).rotate(-angle).map(bar));
  endPainting();
}

void QtImage::drawArrowD(const double midX, const double midY, const double deltaX, const double deltaY, const Color &col, const double size) {
  const double HeadLenFraction = 1.0;
  const double HeadWidthFraction = 1.0;
  const double BarWidthFraction = 0.3;
  const double CircleRadiusFraction = 0.3;

  QPainter *painter = startPainting();
  painter->setRenderHint(QPainter::Antialiasing, _antialiasing);
  painter->setPen(QPen(Qt::NoPen));
  painter->setBrush(QBrush(MyQColor(col)));

  double len = sqrt(deltaX*deltaX + deltaY*deltaY);
  double headLen = HeadLenFraction*size;
  double halfHeadWidth = 0.5*HeadWidthFraction*size;
  double halfBarWidth = 0.5*BarWidthFraction*size;
  double angle = atan2(deltaY, deltaX)/M_PI*180.0;

  if(len<headLen) {
    double circleRadius = CircleRadiusFraction*size;
    painter->drawEllipse(QPointF(midX, yTransformation(midY)), circleRadius, circleRadius);
    endPainting();
    return;
  }

  QPointF firstPoint(-0.5*len, -halfBarWidth);
  QVector<QPointF> v;
  v.append(firstPoint);
  v.append(QPointF(-0.5*len, halfBarWidth));
  v.append(QPointF(0.5*len-headLen, halfBarWidth));
  v.append(QPointF(0.5*len-headLen, halfHeadWidth));
  v.append(QPointF(0.5*len, 0.0));
  v.append(QPointF(0.5*len-headLen, -halfHeadWidth));
  v.append(QPointF(0.5*len-headLen, -halfBarWidth));
  v.append(firstPoint);
  QPolygonF arrow(v);
  painter->drawPolygon(QMatrix().translate(midX, yTransformation(midY)).rotate(-angle).map(arrow));
  endPainting();
}

void QtImage::drawPolygon(const std::vector<Vector2D> &vertices, const Color &col) {
  QPainter *painter = startPainting();
  painter->setRenderHint(QPainter::Antialiasing, _antialiasing);
  painter->setPen(QPen(Qt::NoPen));
  painter->setBrush(QBrush(MyQColor(col)));

  QVector<QPointF> v;
  for(unsigned int i=0; i<vertices.size(); ++i) {
    v.append(QPointF(vertices[i].x(), yTransformation(vertices[i].y())));
  }
  QPolygonF poly(v);
  painter->drawPolygon(poly);
  endPainting();
}  

void QtImage::drawPolygonWithOutline(const std::vector<Vector2D> &vertices, const Color &col, const double width, const Color &outerCol) {
  QPainter *painter = startPainting();
  painter->setRenderHint(QPainter::Antialiasing, _antialiasing);
  painter->setPen(QPen(QBrush(MyQColor(outerCol)), width));
  painter->setBrush(QBrush(MyQColor(col)));

  QVector<QPointF> v;
  for(unsigned int i=0; i<vertices.size(); ++i) {
    v.append(QPointF(vertices[i].x(), yTransformation(vertices[i].y())));
  }
  QPolygonF poly(v);
  painter->drawPolygon(poly);
  endPainting();
}

void QtImage::drawArrowExt(const double midX, const double midY, const double deltaX, const double deltaY, const Color &col, const Color &headCol, const double size) {
  const double HeadLenFraction = 1.0;
  const double HeadWidthFraction = 1.0;
  const double BarWidthFraction = 0.35;
  const double CircleRadiusFraction = 0.3;

  QPainter *painter = startPainting();
  painter->setRenderHint(QPainter::Antialiasing, _antialiasing);
  painter->setPen(QPen(Qt::NoPen));

  double len = sqrt(deltaX*deltaX + deltaY*deltaY);
  double headLen = HeadLenFraction*size;
  double halfHeadWidth = 0.5*HeadWidthFraction*size;
  double halfBarWidth = 0.5*BarWidthFraction*size;
  double angle = atan2(deltaY, deltaX)/M_PI*180.0;

  if(len<headLen) {
    double circleRadius = CircleRadiusFraction*size;
    painter->setBrush(QBrush(MyQColor(headCol)));
    painter->drawEllipse(QPointF(midX, yTransformation(midY)), circleRadius, circleRadius);
    endPainting();
    return;
  }

  painter->setBrush(QBrush(MyQColor(col)));
  QPointF firstPointBar(-0.5*len, -halfBarWidth);
  QVector<QPointF> vBar;
  vBar.append(firstPointBar);
  vBar.append(QPointF(-0.5*len, halfBarWidth));
  vBar.append(QPointF(0.5*len-headLen, halfBarWidth));
  vBar.append(QPointF(0.5*len-headLen, -halfBarWidth));
  vBar.append(firstPointBar);
  QPolygonF bar(vBar);
  painter->drawPolygon(QMatrix().translate(midX, yTransformation(midY)).rotate(-angle).map(bar));

  painter->setBrush(QBrush(MyQColor(headCol)));
  QPointF firstPointHead(0.5*len-headLen, halfHeadWidth);
  QVector<QPointF> vHead;
  vHead.append(firstPointHead);
  vHead.append(QPointF(0.5*len, 0.0));
  vHead.append(QPointF(0.5*len-headLen, -halfHeadWidth));
  vHead.append(firstPointHead);
  QPolygonF head(vHead);
  painter->drawPolygon(QMatrix().translate(midX, yTransformation(midY)).rotate(-angle).map(head));

  endPainting(); 
}

void QtImage::drawArrowExtC(const double midX, const double midY, const double deltaX, const double deltaY, const Color &col, const Color &headCol, const Color &cCol, const double size) {
  const double HeadLenFraction = 1.0;
  const double HeadWidthFraction = 1.0;
  const double BarWidthFraction = 0.35;
  const double CircleRadiusFraction = 0.3;

  QPainter *painter = startPainting();
  painter->setRenderHint(QPainter::Antialiasing, _antialiasing);
  painter->setPen(QPen(Qt::NoPen));

  double len = sqrt(deltaX*deltaX + deltaY*deltaY);
  double headLen = HeadLenFraction*size;
  double halfHeadWidth = 0.5*HeadWidthFraction*size;
  double halfBarWidth = 0.5*BarWidthFraction*size;
  double angle = atan2(deltaY, deltaX)/M_PI*180.0;

  if(len>=headLen) {
    painter->setBrush(QBrush(MyQColor(col)));
    QPointF firstPointBar(-0.5*len, -halfBarWidth);
    QVector<QPointF> vBar;
    vBar.append(firstPointBar);
    vBar.append(QPointF(-0.5*len, halfBarWidth));
    vBar.append(QPointF(0.5*len-headLen, halfBarWidth));
    vBar.append(QPointF(0.5*len-headLen, -halfBarWidth));
    vBar.append(firstPointBar);
    QPolygonF bar(vBar);
    painter->drawPolygon(QMatrix().translate(midX, yTransformation(midY)).rotate(-angle).map(bar));

    painter->setBrush(QBrush(MyQColor(headCol)));
    QPointF firstPointHead(0.5*len-headLen, halfHeadWidth);
    QVector<QPointF> vHead;
    vHead.append(firstPointHead);
    vHead.append(QPointF(0.5*len, 0.0));
    vHead.append(QPointF(0.5*len-headLen, -halfHeadWidth));
    vHead.append(firstPointHead);
    QPolygonF head(vHead);
    painter->drawPolygon(QMatrix().translate(midX, yTransformation(midY)).rotate(-angle).map(head));
  }

  double circleRadius = CircleRadiusFraction*size;
  painter->setBrush(QBrush(MyQColor(cCol)));
  painter->drawEllipse(QPointF(midX, yTransformation(midY)), circleRadius, circleRadius);
  
  endPainting(); 
}

void QtImage::drawText(const double x, const double y, const double pixelSize, const Color &col, const char* text) {
  if(!QApplication::instance()) {
    printf("QtRasterImage::drawText: You should create a QApplication in order to draw text!\n");
    exit(1);
  }
  QPainter *painter = startPainting();
  QFont font(painter->font()); 
  font.setPixelSize(pixelSize);
  painter->setPen(MyQColor(col));
  painter->setFont(font);
  painter->drawText(QPointF(x, yTransformation(y)), text);
  endPainting();
}

void QtImage::drawLineD(const double x1, const double y1, const double x2, const double y2, const Color &col) {
  QPainter *painter = startPainting();
  painter->setRenderHint(QPainter::Antialiasing, _antialiasing);
  painter->setPen(MyQColor(col));
  painter->drawLine(QPointF(x1, yTransformation(y1)), QPointF(x2, yTransformation(y2)));
  endPainting();
}

void QtImage::drawCircle(const double x, const double y, const double radius, const Color &col) {
  QPainter *painter = startPainting();
  painter->setRenderHint(QPainter::Antialiasing, _antialiasing);
  painter->setPen(QPen(Qt::NoPen));
  painter->setBrush(QBrush(MyQColor(col)));
  painter->drawEllipse(QPointF(x, yTransformation(y)), radius, radius);
  endPainting();
}

void QtImage::drawEmptyCircle(const double x, const double y, const double radius, const Color &col, const double thickness) {
  QPainter *painter = startPainting();
  painter->setRenderHint(QPainter::Antialiasing, _antialiasing);
  painter->setPen(QPen(QBrush(MyQColor(col)), thickness));
//   painter->setBrush(QBrush(QColor(r, g, b)));
  painter->drawEllipse(QPointF(x, yTransformation(y)), radius, radius);
  endPainting();
}

void QtImage::drawEllipse(const double x, const double y, const double radius1, const double radius2, const double angle, const Color &col) {
  QPainter *painter = startPainting();
  painter->setRenderHint(QPainter::Antialiasing, _antialiasing);
  painter->setPen(QPen(Qt::NoPen));
  painter->setBrush(QBrush(MyQColor(col)));
  
  painter->translate(x, yTransformation(y));
  painter->rotate(-180.0*angle/M_PI);
  painter->drawEllipse(QPointF(0,0), radius1, radius2);
  
  endPainting();
}
