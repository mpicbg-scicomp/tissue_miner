/* 
 * File:   PixelFrame.cpp
 * Author: mmpi
 * 
 * Created on September 13, 2013, 12:46 PM
 */

#include "PixelFrame.h"

const double Sqrt2 = sqrt(2.0);

// first direction is in positive x direction
const int PixelFrame::NeighborOffsetX[PixelFrame::NumberOfNeighbors] = { 1, 1, 0, -1, -1, -1, 0, 1};
const int PixelFrame::NeighborOffsetY[PixelFrame::NumberOfNeighbors] = { 0, 1, 1, 1, 0, -1, -1, -1};
const double PixelFrame::NNeighborDistance[PixelFrame::NumberOfNeighbors] = { 1, Sqrt2, 1, Sqrt2, 1, Sqrt2, 1, Sqrt2};
