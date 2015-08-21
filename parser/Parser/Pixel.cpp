#include <iostream>

#include "PixelFrame.h"
#include "Pixel.h"

const PixelValue Pixel::BondValue           = 0x00FFFFFF;
const PixelValue Pixel::OutsideCanvasValue  = 0xFF000000;
const PixelValue Pixel::DividingCellValue   = 0x000000FF;

int Pixel::getLastCommonNeighboringBondPixelInCwOrientation(std::vector<PixelValue> &neighbors, const PixelValue v, const int oldDirection) const {
  // clear neighbor set
  neighbors.clear();

  // start from the direction from where I come
  const int directionOffset = PixelFrame::oppositeDirection(oldDirection);

  // get first bond neighbor in cw direction after cell v starting from directionOffset
  int firstBondIndexAfterV=-1;
  PixelValue lastValue = neighbor(directionOffset).data();
  // go around in cw direction
  for(int index=PixelFrame::NumberOfNeighbors-1; index>=0; --index) {
    int curIndex = (directionOffset+index)%PixelFrame::NumberOfNeighbors;
    PixelValue curValue = isNeighborOnCanvas(curIndex)?neighbor(curIndex).data():OutsideCanvasValue;
    if(curValue!=lastValue) {
      // pixel value changes
      if(curValue==BondValue) {
        // change from cell to bond
        if((lastValue==v) && (firstBondIndexAfterV<0)) {
          firstBondIndexAfterV = curIndex;
        }
        neighbors.push_back(lastValue);
      } else if(lastValue!=BondValue) {
        std::cout << "Pixel::getLastCommonNeighboringBondPixelInCwOrientation: Two different cells without bond in between!" << std::endl;
        return -1;        
      }
      lastValue = curValue;
    }
  }
  
  // when I'm here, everything went fine, so far
  // but did we find anything?
  if(firstBondIndexAfterV<0) {
    std::cout << "Pixel::getLastCommonNeighboringBondPixelInCwOrientation: Pixel value not found among neighbors!" << std::endl;
    return -1;
  }
  
  // now, check if we have to take the cw neighbor of firstBondIndexAfterV
  // it is common neighbor of *this and the cell with index v, iff it is not diagonal
  int indexCwNeighbor = (PixelFrame::NumberOfNeighbors-1+firstBondIndexAfterV)%PixelFrame::NumberOfNeighbors; 
  if(!PixelFrame::diagonalNeighbor(indexCwNeighbor) && (neighbor(indexCwNeighbor).data()==BondValue)) {
    return indexCwNeighbor;
  } else {
    // first neighbor in cw direction must be common neighbor and must be a bond pixel (due to consistency check above)
    return firstBondIndexAfterV;
  }
}


void Pixel::addToThickVertex(std::set<Pixel> &thickVertex) {
  thickVertex.insert(*this);
  for(int index=0; index<PixelFrame::NumberOfNeighbors; index+=2) {
    Pixel n1(neighbor(index));
    Pixel n2(neighbor((index+1)%PixelFrame::NumberOfNeighbors));
    Pixel n3(neighbor((index+2)%PixelFrame::NumberOfNeighbors));
    if(!isOnMargin() && (n1.data()==BondValue) && (n2.data()==BondValue) && (n3.data()==BondValue)) {
      if(!thickVertex.count(n1)) {
        n1.addToThickVertex(thickVertex);
      }
      if(!thickVertex.count(n2)) {
        n2.addToThickVertex(thickVertex);
      }
      if(!thickVertex.count(n3)) {
        n3.addToThickVertex(thickVertex);
      }
    }
  }
}
