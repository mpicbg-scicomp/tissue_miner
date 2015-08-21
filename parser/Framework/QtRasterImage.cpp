#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <Qt/QtCore>
#include <Qt/QtGui>

#include "QtRasterImage.h"

QRgb myQRgb(const Color &col) {
  return qRgb(col.r(), col.g(), col.b());
}

QtRasterImage::QtRasterImage() {
  _image = new QImage();
  _valid = false;
  setAlphaChannel(0xFF);
}

QtRasterImage::QtRasterImage(const char filename[]) {
  _image = new QImage();
  _valid = false;
  setAlphaChannel(0xFF);
  load(filename);
}

QtRasterImage::QtRasterImage(const int width, const int height) {
  _image = new QImage(QSize(width, height), QImage::Format_ARGB32);
  setAlphaChannel(0xFF);
  _valid = true;
}

QtRasterImage::~QtRasterImage() {
  delete _image;
}

double QtRasterImage::yTransformation(const double y) const {
  return _image->height()-1-y;
}

bool QtRasterImage::load(const char filename[]) {
  if(!_image->load(filename)) {
    printf("QtRasterImage::load: Could not load image %s!\n", filename);
    _valid = false;
    return false;
  }
  *_image = _image->convertToFormat(QImage::Format_ARGB32);
  _valid = true;
  return true;
}

bool QtRasterImage::save(const char filename[]) const {
  if(!_image->save(filename)) {
    printf("QtRasterImage::save: Could not save image %s!\n", filename);
    QImageWriter writer(filename, "PNG");
    if(!writer.write(*_image)) {
      printf("QtRasterImage::save: Could not save image %s via QImageWriter class: %s --- %s\n", filename, writer.errorString().toLatin1().data(), writer.device()->errorString().toLatin1().data());
    }
    return false;
  }
  printf("Saved %s.\n", filename);
  return true;
}

void QtRasterImage::reset(const Color &col) {
  _image->fill(myQRgb(col));
}

void QtRasterImage::invert() {
  int width=_image->width(), height=_image->height();
  QRgb* scanLine;
  for(int j=0; j<height; ++j) {
    scanLine = (QRgb*)_image->scanLine(height-j-1);
    for(int i=0; i<width; ++i) {
      scanLine[i] = qRgba(255-qRed(scanLine[i]), 255-qGreen(scanLine[i]), 255-qBlue(scanLine[i]), qAlpha(scanLine[i]));
    }
  }
}

void QtRasterImage::setRedChannel(const unsigned char r) {
  int width=_image->width(), height=_image->height();
  QRgb* scanLine;
  for(int j=0; j<height; ++j) {
    scanLine = (QRgb*)_image->scanLine(height-j-1);
    for(int i=0; i<width; ++i) {
      scanLine[i] = qRgba(r, qGreen(scanLine[i]), qBlue(scanLine[i]), qAlpha(scanLine[i]));
    }
  }
}

void QtRasterImage::setGreenChannel(const unsigned char g) {
  int width=_image->width(), height=_image->height();
  QRgb* scanLine;
  for(int j=0; j<height; ++j) {
    scanLine = (QRgb*)_image->scanLine(height-j-1);
    for(int i=0; i<width; ++i) {
      scanLine[i] = qRgba(qRed(scanLine[i]), g, qBlue(scanLine[i]), qAlpha(scanLine[i]));
    }
  }
}

void QtRasterImage::setBlueChannel(const unsigned char b) {
  int width=_image->width(), height=_image->height();
  QRgb* scanLine;
  for(int j=0; j<height; ++j) {
    scanLine = (QRgb*)_image->scanLine(height-j-1);
    for(int i=0; i<width; ++i) {
      scanLine[i] = qRgba(qRed(scanLine[i]), qGreen(scanLine[i]), b, qAlpha(scanLine[i]));
    }
  }
}

void QtRasterImage::setAlphaChannel(const unsigned char a) {
  int width=_image->width(), height=_image->height();
  QRgb* scanLine;
  for(int j=0; j<height; ++j) {
    scanLine = (QRgb*)_image->scanLine(height-j-1);
    for(int i=0; i<width; ++i) {
      scanLine[i] = qRgba(qRed(scanLine[i]), qGreen(scanLine[i]), qBlue(scanLine[i]), a);
    }
  }
}

void QtRasterImage::copyToRedChannel(const unsigned char *const data) {
  int width=_image->width(), height=_image->height();
  QRgb* scanLine;
  for(int j=0; j<height; ++j) {
    scanLine = (QRgb*)_image->scanLine(height-j-1);
    for(int i=0; i<width; ++i) {
//      scanLine[i] = qRbg(qRed(scanLine[i]), qGreen(scanLine[i]), qBlue(scanLine[i]));
      scanLine[i] = qRgba(data[i*height+j], qGreen(scanLine[i]), qBlue(scanLine[i]), qAlpha(scanLine[i]));
    }
  }
}

void QtRasterImage::copyToGreenChannel(const unsigned char *const data) {
  int width=_image->width(), height=_image->height();
  QRgb* scanLine;
  for(int j=0; j<height; ++j) {
    scanLine = (QRgb*)_image->scanLine(height-j-1);
    for(int i=0; i<width; ++i) {
//      scanLine[i] = qRbg(qRed(scanLine[i]), qGreen(scanLine[i]), qBlue(scanLine[i]));
      scanLine[i] = qRgba(qRed(scanLine[i]), data[i*height+j], qBlue(scanLine[i]), qAlpha(scanLine[i]));
    }
  }
}

void QtRasterImage::copyToBlueChannel(const unsigned char *const data) {
  int width=_image->width(), height=_image->height();
  QRgb* scanLine;
  for(int j=0; j<height; ++j) {
    scanLine = (QRgb*)_image->scanLine(height-j-1);
    for(int i=0; i<width; ++i) {
//      scanLine[i] = qRbg(qRed(scanLine[i]), qGreen(scanLine[i]), qBlue(scanLine[i]));
      scanLine[i] = qRgba(qRed(scanLine[i]), qGreen(scanLine[i]), data[i*height+j], qAlpha(scanLine[i]));
    }
  }
}

void QtRasterImage::copyToAlphaChannel(const unsigned char *const data) {
  int width=_image->width(), height=_image->height();
  QRgb* scanLine;
  for(int j=0; j<height; ++j) {
    scanLine = (QRgb*)_image->scanLine(height-j-1);
    for(int i=0; i<width; ++i) {
//      scanLine[i] = qRbg(qRed(scanLine[i]), qGreen(scanLine[i]), qBlue(scanLine[i]));
      scanLine[i] = qRgba(qRed(scanLine[i]), qGreen(scanLine[i]), qBlue(scanLine[i]), data[i*height+j]);
    }
  }
}

unsigned char* QtRasterImage::redChannel() const {
  int width=_image->width(), height=_image->height();
  unsigned char *data = new unsigned char[width*height];
  QRgb* scanLine;
  for(int j=0; j<height; ++j) {
    scanLine = (QRgb*)_image->scanLine(height-j-1);
    for(int i=0; i<width; ++i) {
      data[j+i*height] = qRed(scanLine[i]);
    }
  }
  return data;
}

unsigned char* QtRasterImage::greenChannel() const {
  int width=_image->width(), height=_image->height();
  unsigned char *data = new unsigned char[width*height];
  QRgb* scanLine;
  for(int j=0; j<height; ++j) {
    scanLine = (QRgb*)_image->scanLine(height-j-1);
    for(int i=0; i<width; ++i) {
      data[j+i*height] = qGreen(scanLine[i]);
    }
  }
  return data;
}

unsigned char* QtRasterImage::blueChannel() const {
  int width=_image->width(), height=_image->height();
  unsigned char *data = new unsigned char[width*height];
  QRgb* scanLine;
  for(int j=0; j<height; ++j) {
    scanLine = (QRgb*)_image->scanLine(height-j-1);
    for(int i=0; i<width; ++i) {
      data[j+i*height] = qBlue(scanLine[i]);
    }
  }
  return data;
}

unsigned char* QtRasterImage::alphaChannel() const {
  int width=_image->width(), height=_image->height();
  unsigned char *data = new unsigned char[width*height];
  QRgb* scanLine;
  for(int j=0; j<height; ++j) {
    scanLine = (QRgb*)_image->scanLine(height-j-1);
    for(int i=0; i<width; ++i) {
      data[j+i*height] = qAlpha(scanLine[i]);
    }
  }
  return data;
}

unsigned int* QtRasterImage::rgbaData() const {
  int width=_image->width(), height=_image->height();
  unsigned int *data = new unsigned int[width*height];
  QRgb* scanLine;
  for(int j=0; j<height; ++j) {
    scanLine = (QRgb*)_image->scanLine(height-j-1);
    for(int i=0; i<width; ++i) {
      data[j+i*height] = (unsigned int)scanLine[i];
    }
  }
  return data;
}

unsigned int* QtRasterImage::rgbData() const {
  int width=_image->width(), height=_image->height();
  unsigned int *data = new unsigned int[width*height];
  QRgb* scanLine;
  for(int j=0; j<height; ++j) {
    scanLine = (QRgb*)_image->scanLine(height-j-1);
    for(int i=0; i<width; ++i) {
      data[j+i*height] = 0x00FFFFFF & (unsigned int)scanLine[i];
    }
  }
  return data;
}

void QtRasterImage::drawPixel(const int i, const int j, const Color &col) {
  _image->setPixel(i, _image->height()-j-1, myQRgb(col));
}

char QtRasterImage::pixel(const int i, const int j, const ColorChannel c) const {
  if(c==RedChannel) {
    return qRed(_image->pixel(i, _image->height()-j-1));
  } else if(c==GreenChannel) {
    return qGreen(_image->pixel(i, _image->height()-j-1));
  } else if(c==BlueChannel) {
    return qBlue(_image->pixel(i, _image->height()-j-1));
  } else if(c==AlphaChannel) {
    return qAlpha(_image->pixel(i, _image->height()-j-1));
  }
  return 0;
}

QtRasterImage QtRasterImage::scaled(int width, int height) {
  QtRasterImage img;
  img._image = new QImage(_image->scaled(width, height, Qt::KeepAspectRatio, Qt::SmoothTransformation));
  img._valid = false;
  return img;
}
