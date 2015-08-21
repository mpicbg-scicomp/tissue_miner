#include <stdio.h>
#include <math.h>

#include "AbstractCanvas.h"


// adopted from Doug
void AbstractCanvas::drawArrowD(const double midX, const double midY, const double deltaX, const double deltaY, const Color &col, const double size) {
  double length = sqrt( deltaX*deltaX + deltaY*deltaY );
  if(length>=0.5*size) {
    double dXdS = deltaX/length;
    double dYdS = deltaY/length;
  
    double headX = midX+0.5*deltaX;
    double headY = midY+0.5*deltaY;
    int cwX  = (int) (headX-size*dXdS+size*dYdS);
    int cwY  = (int) (headY-size*dYdS-size*dXdS);
    int ccwX = (int) (headX-size*dXdS-size*dYdS);
    int ccwY = (int) (headY-size*dYdS+size*dXdS);
    drawLineD(headX, headY, cwX, cwY, col);
    drawLineD(headX, headY, ccwX, ccwY, col);
  
    double tailX = midX-0.5*deltaX;
    double tailY = midY-0.5*deltaY;
    drawLineD(headX, headY, tailX, tailY, col);
  } else {
    drawLineD(midX, midY, midX, midY, col);
  }
}

void AbstractCanvas::drawBarD(const double x1, const double y1, const double x2, const double y2, const Color &col, const double size) {
  drawLineD(x1, y1, x2, y2, col);
}

