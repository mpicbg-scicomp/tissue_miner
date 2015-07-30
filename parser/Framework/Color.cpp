#include <stdio.h>
#include <math.h>

#include "Color.h"

const Color Black(0x00, 0x00, 0x00);
const Color White(0xFF, 0xFF, 0xFF);
const Color Red(0xFF, 0x00, 0x00);
const Color Green(0x00, 0xFF, 0x00);
const Color Blue(0x00, 0x00, 0xFF);
const Color Cyan(0x00, 0xFF, 0xFF);
const Color Magenta(0xFF, 0x00, 0xFF);
const Color Yellow(0xFF, 0xFF, 0x00);

const Color Orange(0xFF, 0x80, 0x00);
const Color Lightblue(0x40, 0x40, 0xFF);


bool ColorBar::saveAsLut(const char *fileName) const {
  FILE *f = fopen(fileName, "w");
  if(!f) {
    return false;
  }
  for(int i=0; i<256; ++i) {
    Color col(color(-1.0+2.0*i/255.0));
    fprintf(f, "%d %d %d\n", (int)col.r(), (int)col.g(), (int)col.b());
  }
  fclose(f);
  return true;
};

Color RainbowColorBar::color(const double v) const {
  double v2 = (v<1)?( (v>-1)?v:-1 ):1;
  unsigned char r = (v2>0)?255*v2:0;
  unsigned char g = (1.0-fabs(v2))*255;
  unsigned char b = (v2<0)?-255*v2:0;
  return Color(r, g, b);  
};

LutColorBar::LutColorBar(const char *fileName) : _valid(false) {
  FILE *f = fopen(fileName, "r");
  if(!f) {
    return;
  }
  for(int i=0; i<256; ++i) {
    if(feof(f)) {
      return;
    }
    int r, g, b;
    if(3!=fscanf(f, " %d %d %d", &r, &g, &b)) {
      return;
    }
    _r[i] = (unsigned char)r;
    _g[i] = (unsigned char)g;
    _b[i] = (unsigned char)b;
  }
  fclose(f);
  _valid = true;
}

Color LutColorBar::color(const double v) const {
  int i = int((v+1.0)/2.0*255.0);
  return Color(_r[i], _g[i], _b[i]);
}

const RainbowColorBar Rainbow;

