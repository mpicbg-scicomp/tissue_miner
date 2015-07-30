/* 
 * File:   Cell.h
 * Author: mmpi
 *
 * Created on September 13, 2013, 11:00 AM
 */

#ifndef CELL_H
#define	CELL_H

#include <vector>

#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif

#include "DirectedBond.h"

/** cell index type */
typedef unsigned int CellIndex;

/** contains cell properties and references to bonds */
struct Cell {
  Cell(const CellIndex id_) : id(id_), duringTransitionBefore(UnclassifiedBefore), mother(0), sister(0), duringTransitionAfter(UnclassifiedAfter), daughter(0), area(0), r(0,0), elongation(0,0), deformedCell(NULL) {};
  Cell(const Cell &other) : id(other.id), duringTransitionBefore(other.duringTransitionBefore), mother(other.mother), sister(other.sister), duringTransitionAfter(other.duringTransitionAfter), daughter(other.daughter), area(other.area), r(other.r), deformedCell(NULL) {};

#ifdef USE_NETCDF
  /** for saving the data contained in all cell objects */
  class DataFileIO {
  public:
    DataFileIO() : _numberOfFields(0) {}
    ~DataFileIO();
    bool load(TissueState &s, const NcFile &dataFile);
    bool save(NcFile &dataFile, const TissueState &s);
  private:
    void allocateBuffers(int numberOfFields);
    int _numberOfFields;
    NcVar *_idVar, *_areaVar, *_xVar, *_yVar, *_duringTransitionBeforeVar, *_motherVar, *_sisterVar, *_duringTransitionAfterVar, *_daughterVar, *_e1Var, *_e2Var, *_pR1Var, *_pR2Var, *_iRVar, *_pG1Var, *_pG2Var, *_iGVar, *_pB1Var, *_pB2Var, *_iBVar;
    int *_id, *_duringTransitionBefore, *_mother, *_sister, *_duringTransitionAfter, *_daughter;
    double *_area, *_x, *_y, *_e1, *_e2, *_pR1, *_pR2, *_iR, *_pG1, *_pG2, *_iG, *_pB1, *_pB2, *_iB;
  };
#endif  
  
  /** This is the cell id (encoded in the color, format 0xAARRGGBB). */
  const CellIndex id;
  static const CellIndex VoidCellId = 0;
  
  // during transition before
  enum {
    UnclassifiedBefore=0, Stayed, Divided, MovedIntoMask, SegmentationErrorAppearance
  } duringTransitionBefore; /*< what happened to this cell before this frame? */
  /** mother cell id */
  CellIndex mother;
  /** sister cell id */
  CellIndex sister;
  
  // during transition after
  enum {
    UnclassifiedAfter=0, Stays, Divides, Apoptosis, MovesOutOfMask, SegmentationErrorDisappearance
  } duringTransitionAfter; /*< what happens to this cell after this frame? */
  /** daughter with different cell id */
  CellIndex daughter;
  
  /** These are pointers to all bonds of this cell sorted in ccw order. The directed bonds point into ccw direction. */
  std::vector<DirectedBond *> bonds;
  
  /** cell area using to irregular pixel outline */
  double area;
  
  /** cell center using to irregular pixel outline; area center of cell */
  Vector2D r;
  
  /** cell elongation according to cell paper definition */
  Nematic2D elongation;
  
  /** cell polarity according to cell paper definition */
  Nematic2D polarityR, polarityG, polarityB;
  double intIntensityR, intIntensityG, intIntensityB;
  
  /** this is a pointer to the (moved) copy of this cell */
  Cell *deformedCell;
};

#endif	/* CELL_H */

