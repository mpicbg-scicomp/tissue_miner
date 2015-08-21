/* 
 * File:   DirectedBond.h
 * Author: mmpi
 *
 * Created on September 13, 2013, 11:00 AM
 */

#ifndef DIRECTEDBOND_H
#define	DIRECTEDBOND_H

#include <map>

class Cell;
class Vertex;
class TissueState;

/**
 * A directed bond is referenced by rightVertex, cell is left seen from rightVertex.
 * A directed bond is referenced by cell, rightVertex is right seen from the cell.
 */
struct DirectedBond {
  DirectedBond() : conjBond(NULL), cell(NULL), rightVertex(NULL), leftVertex(NULL), length(0.0) {}

  /** when valid state, if this is NULL, then, this is a margin bond */
  DirectedBond *conjBond;
  /** when valid state, this is never NULL */
  Cell *cell;
  /** when valid state, this is never NULL */
  Vertex *rightVertex, *leftVertex;  
  
  /** length of (curved) bond in pixels, superposition of 1 and sqrt(2) */
  double length;
};


#endif	/* DIRECTEDBOND_H */

