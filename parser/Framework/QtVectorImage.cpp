#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <Qt/QtCore>
#include <Qt/QtGui>

#include "Transformation.h"
#include "QtRasterImage.h"
#include "QtVectorImage.h"


void QtVectorImage::init(const int width, const int height) {
  _printer = new QPrinter();
  _printer->setPaperSize(QSizeF(width, height), QPrinter::DevicePixel);
  _printer->setFullPage(true);
  _printer->setColorMode(QPrinter::Color);
}

QtVectorImage::QtVectorImage(const int width, const int height) : _valid(false) {
  init(width, height);
}

QtVectorImage::QtVectorImage(const QtRasterImage &img) : _valid(false) {
  init(img.width(), img.height());
}

QtVectorImage::~QtVectorImage() {
  delete _printer;
}

double QtVectorImage::yTransformation(const double y) const {
  return _printer->height()-y;
}

void QtVectorImage::reset(const Color &col) {
  QPainter *painter = startPainting();
  painter->setPen(QPen(Qt::NoPen)); //QColor(r, g, b)
  painter->setBrush(QBrush(MyQColor(col)));
  painter->drawRect(0.0, 0.0, width(), height());
  endPainting();
}

void QtVectorImage::drawImage(const QtRasterImage &img, const double x, const double y) {
  QPainter *painter = startPainting();
  painter->drawImage(x, height()-y-img.height(), *img._image);
  endPainting();
}

bool QtVectorImage::start(const char filename[]) {
  _printer->setOutputFileName(filename);
  _painter = new QPainter();
  _valid = _painter->begin(_printer);
  if(!_valid) {
    _painter->end();
    delete _painter;
    printf("Could not open %s for writing!.\n", _printer->outputFileName().toAscii().constData());
  }
  return _valid;
}

void QtVectorImage::done() {
  if(_valid) {
    _painter->end();
    delete _painter;
    printf("Saved %s.\n", _printer->outputFileName().toAscii().constData());
  }
}

