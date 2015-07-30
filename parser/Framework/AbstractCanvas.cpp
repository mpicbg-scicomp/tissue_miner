#include <stdio.h>
#include <math.h>

#include "AbstractCanvas.h"
#include "AbstractField.h"


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


void AbstractCanvas::drawNematic2DArray(const IsotropicTransformation &trafo, const AbstractNematic2DArray &a, const bool normalized, const Color &Col, const double RelativeLength, const double RelWidth) {
  const double NematicLength = RelativeLength*a.grid().minBoxSize()*fabs(trafo.scale());
  const double BarSize = RelWidth*a.grid().minBoxSize()*fabs(trafo.scale());
  for(int i=0; i<a.grid().numBoxesX(); ++i) {
    for(int j=0; j<a.grid().numBoxesY(); ++j) {
      if(a.valid(i,j)) {
        Nematic2D n(a(i,j));
        n = nematicFromAngleAndNorm(trafo.mapAngle(n.angle()), normalized?1.0:n.norm());
        drawBar(trafo.map(a.grid().boxMid(i,j)), NematicLength*n, Col, BarSize);
      }
    }
  }  
}
