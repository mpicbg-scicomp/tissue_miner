#include <math.h>
#include "AbstractRasterCanvas.h"

// adopted from Doug
void AbstractRasterCanvas::drawLineD(const double x1, const double y1, const double x2, const double y2, const Color &col) {
  double deltaX = x2-x1;
  double deltaY = y2-y1;
  double length = sqrt( deltaX*deltaX + deltaY*deltaY );
  
  double dXdS = deltaX/length;
  double dYdS = deltaY/length;
  for(int s=0; s<=length; s++){
    int i = (int)(x1+dXdS*s);
    int j = (int)(y1+dYdS*s);
    if(i>=0 && j>=0 && i<width() && j<height())
      drawPixel(i, j, col);
  }
}

unsigned char* AbstractRasterCanvas::channel(ColorChannel c) const {
  switch(c) {
    case RedChannel:
      return redChannel();
    case GreenChannel:
      return greenChannel();
    case BlueChannel:
      return blueChannel();
    case AlphaChannel:
      return alphaChannel();
  }
  return NULL;
}
