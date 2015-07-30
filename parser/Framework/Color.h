#ifndef COLOR_H
#define COLOR_H

class Color {
public:
  Color(const unsigned char r=0x00, const unsigned char g=0x00, const unsigned char b=0x00, const unsigned char a=0xFF) : _r(r), _g(g), _b(b), _a(a) {};
  Color(const Color &other) : _r(other._r), _g(other._g), _b(other._b), _a(other._a) {};
  unsigned char r() const { return _r; } 
  unsigned char g() const { return _g; }
  unsigned char b() const { return _b; }
  unsigned char a() const { return _a; }
  unsigned long rbgToLong() const { return (((((unsigned long)_r) << 8) + _g) << 8) + _b; }
  
  bool operator==(Color other) {
    return (_r==other._r) && (_g==other._g) && (_b==other._b) && (_a==other._a);
  }
  
  Color withAlpha(const unsigned int a) const {
    return Color(_r, _g, _b, a);
  }
  
private:  
  unsigned char _r, _g, _b, _a;
};


enum ColorChannel {
  NoChannel=-1, RedChannel, GreenChannel, BlueChannel, AlphaChannel
};

extern const Color Black;
extern const Color White;
extern const Color Red;
extern const Color Green;
extern const Color Blue;
extern const Color Cyan;
extern const Color Magenta;
extern const Color Yellow;

extern const Color Orange;
extern const Color Lightblue;

// --- ColorBar ---

class ColorBar {
public:
  virtual Color color(const double v) const = 0; // v=-1...+1
  bool saveAsLut(const char *fileName) const;
};

class RainbowColorBar : public ColorBar {
public:
  virtual Color color(const double v) const;
};

class LutColorBar : public ColorBar {
public:
  LutColorBar(const char *fileName);
  bool valid() const { return _valid; }
  virtual Color color(const double v) const;
private:
  bool _valid;
  unsigned char _r[256], _g[256], _b[256];
};

extern const RainbowColorBar Rainbow;

#endif

