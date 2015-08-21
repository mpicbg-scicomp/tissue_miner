/* 
 * File:   TissueState.h
 * Author: mmpi
 *
 * Created on September 13, 2013, 11:09 AM
 */

#ifndef TISSUESTATE_H
#define	TISSUESTATE_H

#include <stdlib.h>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <fstream>

#include "PixelFrame.h"
#include "Pixel.h"
#include "Vertex.h"
#include "DirectedBond.h"
#include "Cell.h"

using namespace std;

class TissueState {
public:
  // ********* creation **********
   
  /** creates an empty state at time */
  TissueState(const int frameNumber, const double time) : _frameNumber(frameNumber), _time(time) {}
  ~TissueState();
 
  /**
   * parses tracked_cells image file
   * @param FileName image file name
   * @return true if successful
   */
  bool parseFromTrackedCellsFile(const std::string FileName, const string OriginalFileName);
  
  /** checks leftVertex<->conjBond consistency, bond sorting in vertices, and bond sorting in cells */
  bool checkTopologicalConsistency() const;

  /** removes margin cells */
  void removeMarginCells();
  
  /** exports the whole topology data etc. to text files for database */
  bool exportToDbTables(const int ImageHeight, ofstream& framesFile, ofstream &verticesFile, ofstream &cellsFile, ofstream &ignoredCellsFile, ofstream &undirectedBondsFile, ofstream &directedBondsFile, unsigned long &lastVid, unsigned long &lastDbid, unsigned long &lastUbid) const;

  // ********* element access **********
  
  /** frame number of tracked cells file */
  int frameNumber() const { return _frameNumber; }
  /** time of state */
  double time() const { return _time; }  
  
  /** number of cells in tissue state */
  int numberOfCells() const { return _cells.size(); }
  /** iterator type for cells of a tissue state */
  typedef map<CellIndex, Cell*>::iterator CellIterator;
  /** iterator type for cells of a tissue state */
  typedef map<CellIndex, Cell*>::const_iterator CellConstIterator;
  /** iterator to first cell of tissue */
  CellIterator beginCellIterator() { return _cells.begin(); }
  /** iterator to first cell of tissue */
  CellConstIterator beginCellIterator() const { return _cells.begin(); }
  /** end cell iterator for tissue */
  CellIterator endCellIterator() { return _cells.end(); }
  /** end cell iterator for tissue */
  CellConstIterator endCellIterator() const { return _cells.end(); }
  /** cell for iterator */
  Cell *cell(const CellConstIterator &it) const { return it->second; }
  /** cell for id */
  Cell *cell(const CellIndex id) const { return _cells.at(id); }
  /** check if id is among cell considered */
  bool contains(const CellIndex id) const { return _cells.count(id); }
  /** check if id is among ignored cells on the image margin */
  bool isOnImageMargin(const CellIndex id) const { return _ignoredCells.count(id); }
  /** cell for id, check if id is there */
  Cell *cellCheck(const CellIndex id) const { return contains(id)?_cells.at(id):NULL; }
  /** creates cell with id */
  Cell *newCell(const CellIndex id) {
    if(contains(id)) {
      return NULL;
    } else {
      return _cells[id] = new Cell(id);
    }
  }
  /** entirely removes cell from network */
  void removeCell(CellIndex id);
  
  /** number of bonds in tissue state */
  int numberOfBonds() const { return _bonds.size(); }
  /** iterator type for bonds of a tissue state */
  typedef set<DirectedBond *>::iterator BondIterator;
  /** iterator type for bonds of a tissue state */
  typedef set<DirectedBond *>::const_iterator BondConstIterator;
  /** iterator to first bond of tissue */
  BondIterator beginBondIterator() { return _bonds.begin(); }
  /** iterator to first bond of tissue */
  BondConstIterator beginBondIterator() const { return _bonds.begin(); }
  /** end bond iterator for tissue */
  BondIterator endBondIterator() { return _bonds.end(); }
  /** end bond iterator for tissue */
  BondConstIterator endBondIterator() const { return _bonds.end(); }
  /** bond for iterator */
  DirectedBond *bond(const BondConstIterator &it) const { return *it; }
  /** creates bond */
  DirectedBond *newBond() {
    DirectedBond *b = new DirectedBond;
    _bonds.insert(b);
    return b;
  }
  
  /** number of vertices in tissue state */
  int numberOfVertices() const { return _vertices.size(); }
  /** vertex with index */
  Vertex *vertex(const int index) const { return _vertices[index]; }
  /** creates vertex */
  Vertex *newVertex(const Vector2D &r) {
    Vertex *v = new Vertex(r);
    _vertices.push_back(v);
    return v;
  }

private:
  // ********* properties **********
  /** frame number of tracked cells file */
  int _frameNumber;
  /** time point of state */
  double _time;
  /** holds all vertices of the tissue state */
  vector<Vertex *> _vertices;
  /** holds all bonds of the tissue state */
  set<DirectedBond *> _bonds;
  /** holds all cells of the tissue state */
  map<CellIndex, Cell*> _cells;

  // ********* needed for parsing images *********
  set<CellIndex> _ignoredCells;
  map<Pixel, Vertex*> _vertexMap;
  map<pair<Pixel, int>, DirectedBond*> _directedBondMap;
  /**
   * Creates cell, and goes around cell perimeter of cell in ccw direction. It thereby creates vertices, bonds and the necessary connections.
   * @param id cell id
   * @param firstP boundary pixel
   * @return true if everything is consistent
   */
  bool addCell(const CellIndex id, Pixel firstP, const PixelFrame &OriginalFrame);
          
  // ********* methods **********
  /** entirely removes bond from network; Does not remove bond reference from cell! */
  void removeBondWithoutCell(DirectedBond *b);
  /** empties state */
  void cleanUp();
};

#endif	/* TISSUESTATE_H */

