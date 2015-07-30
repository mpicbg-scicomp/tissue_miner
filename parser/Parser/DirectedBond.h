/* 
 * File:   DirectedBond.h
 * Author: mmpi
 *
 * Created on September 13, 2013, 11:00 AM
 */

#ifndef DIRECTEDBOND_H
#define	DIRECTEDBOND_H

#include <map>

#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif

#include "Triangle.h"

class Cell;
class Vertex;
class TissueState;

/**
 * A directed bond is referenced by rightVertex, cell is left seen from rightVertex.
 * A directed bond is referenced by cell, rightVertex is right seen from the cell.
 */
struct DirectedBond {
  DirectedBond() : conjBond(NULL), cell(NULL), rightVertex(NULL), leftVertex(NULL), deformedBond(NULL), bondTriangle(NULL), length(0.0) {}
  void createTriangle();
  ~DirectedBond();

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
    NcVar *_conjBondVar, *_cellVar, *_rightVertexVar, *_leftVertexVar, *_leftBondSeenFromVertexVar, *_leftBondSeenFromCellVar;
    int *_conjBond, *_cell, *_rightVertex, *_leftVertex, *_leftBondSeenFromVertex, *_leftBondSeenFromCell;
  };
#endif  
  
  /** when valid state, if this is NULL, then, this is a margin bond */
  DirectedBond *conjBond;
  /** when valid state, this is never NULL */
  Cell *cell;
  /** when valid state, this is never NULL */
  Vertex *rightVertex, *leftVertex;  
  
  /** deformed bond */
  DirectedBond *deformedBond;
  
  /** bond triangle */
  Triangle *bondTriangle;

  /** length of (curved) bond in pixels, superposition of 1 and sqrt(2) */
  double length;
};


#endif	/* DIRECTEDBOND_H */

