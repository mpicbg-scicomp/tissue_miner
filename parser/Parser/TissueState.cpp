#include <set>
#include <stack>
#include <vector>
#include <algorithm>
#include <fstream>

#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif

#include "TextProgressBar.h"
#include "QtRasterImage.h"

#include "TissueState.h"
#include "Pixel.h"
#include "LogFile.h"

TissueState::~TissueState() {
  cleanUp();  
}

void TissueState::cleanUp() {
  for(vector<Vertex*>::iterator it=_vertices.begin(); it!=_vertices.end(); ++it) {
    delete *it;
  }
  for(BondIterator it=beginBondIterator(); it!=endBondIterator(); ++it) {
    delete bond(it);
  }
  for(CellIterator it=beginCellIterator(); it!=endCellIterator(); ++it) {
    delete cell(it);    
  }
  _vertices.clear();
  _bonds.clear();
  _cells.clear();
}


void TissueState::removeBondWithoutCell(DirectedBond *b) {
  if(!_bonds.count(b)) {
    std::cout << "TissueState::removeBondWithoutCell: bond not found!" << std::endl;
    throw std::exception();
  }
  // -- bond was found --
  // remove from conj bond
  if(b->conjBond) {
    b->conjBond->conjBond = NULL;
  }
  // remove from vertex
  Vertex *v = b->rightVertex;
  if(!v) {
    std::cout << "TissueState::removeBondWithoutCell: has no rightVertex!" << std::endl;
    throw std::exception();
  }
  std::vector<DirectedBond *>::iterator bondInVertexIt = std::find(v->bonds.begin(), v->bonds.end(), b);
  if(bondInVertexIt==v->bonds.end()) {
    std::cout << "TissueState::removeBondWithoutCell: bond not found within rightVertex!" << std::endl;
    throw std::exception();
  }
  v->bonds.erase(bondInVertexIt);
  // check vertex size
  if(v->bonds.size()==0) {
    std::vector<Vertex*>::iterator vertexIt = std::find(_vertices.begin(), _vertices.end(), v);
    if(vertexIt==_vertices.end()) {
      std::cout << "TissueState::removeBondWithoutCell: vertex not found!" << std::endl;
      throw std::exception();
    }
    delete v;
    _vertices.erase(vertexIt);
  }
  // delete bond   
  delete b;
  _bonds.erase(b);
}

void TissueState::removeBond(DirectedBond *b) {
  if(!_bonds.count(b)) {
    std::cout << "TissueState::removeBond: bond not found!" << std::endl;
    throw std::exception();
  }
  // -- bond was found --
  // remove from cell
  Cell *c = b->cell;
  if(!c) {
    std::cout << "TissueState::removeBond: has no cell!" << std::endl;
    throw std::exception();
  }
  std::vector<DirectedBond *>::iterator bondInCellIt = std::find(c->bonds.begin(), c->bonds.end(), b);
  if(bondInCellIt==c->bonds.end()) {
    std::cout << "TissueState::removeBond: bond not found within cell!" << std::endl;
    throw std::exception();
  }
  c->bonds.erase(bondInCellIt);
  // check cell size
  if(c->bonds.size()==0) {
    removeCell(c->id);
    std::cout << "TissueState::removeBond: removed cell!" << std::endl;
  }
  
  removeBondWithoutCell(b);
}

void TissueState::removeCell(CellIndex id) {
  Cell *c = cell(id);
  for(unsigned int j=0; j<c->bonds.size(); ++j) {
    removeBondWithoutCell(c->bonds[j]);
  }
  delete c;
  _cells.erase(id);
}

bool TissueState::parseFromTrackedCellsFile(const string FileName, const string OriginalFileName) {
  cleanUp();
  
  // load images
  QtRasterImage image(FileName.c_str());
  if(!image.valid()) {
    return false;
  }
  QtRasterImage original(OriginalFileName.c_str());
  if(!original.valid()) {
    return false;
  }
  // prepare pixel frames
  PixelFrame frame(image);
  PixelFrame originalFrame(original);
  
  TextProgressBar bar;
  
  // loop through all pixels of the image
  Pixel p(frame);
  PixelValue lastValue=p.data();
  while(true) {
    Pixel lastP(p);
    ++p;
    // done?
    if(!p.inCanvas()) { break; }
    // get pixel value
    PixelValue curValue=p.data();
    // value changes?
    if(lastValue!=curValue) {
      // was bond before?
      if(lastValue==Pixel::BondValue) {
        // -> curValue corresponds to a cell id and lastP sits on a boundary pixel
        CellIndex id = curValue;
        // cell not yet created and not among the ignored cells?
        if((_cells.count(id)==0) && (_ignoredCells.count(id)==0)) {
          if(!addCell(id, lastP, originalFrame)) {
            return false;
          }
        }
      }
      lastValue = curValue;
    }
    bar.update(p.fractionOfCanvas());
  }
  bar.done(true);
  
  // clear unneeded data
  _vertexMap.clear();
  _directedBondMap.clear();
  
  // remove directed bonds without cell (created in addCell as old conjugated bonds; this was needed for the sorting within vertices)
  std::stack<DirectedBond*> toBeRemoved;
  for(TissueState::BondIterator it=beginBondIterator(); it!=endBondIterator(); ++it) {
    if(!bond(it)->cell) {
      toBeRemoved.push(bond(it));
    }
  }
  while(!toBeRemoved.empty()) {
    removeBondWithoutCell(toBeRemoved.top());
    toBeRemoved.pop();
  }

  
  return true;
}

bool compareCyclicLists(const std::vector<CellIndex> &list1, const std::vector<CellIndex> &list2) {
  if(list1.size()!=list2.size()) {
    return false;
  }
  // look for element in list 2 that corresponds to first element in list1
  int found = -1;
  for(unsigned int i=0; i<list2.size(); ++i) {
    if(list2[i]==list1[0]) {
      found = i;
      break;
    }
  }
  if(found<0) {
    return false;
  }
  for(unsigned int i=0; i<list1.size(); ++i) {
    if(list1[i]!=list2[(i+found)%list2.size()]) {
      return false;
    }
  }
  return true;
}

bool TissueState::addCell(const CellIndex id, const Pixel firstP, const PixelFrame &OriginalFrame) {
  // create cell, already now for computation of geometric quantities...
  Cell *c = new Cell(id);

  // **** first, check if cell has to be ignored (i.e. at least one bond pixel is a margin pixel) ****
  // go around cell perimeter in ccw direction
  Pixel p(firstP);
  Vector2D lastR(p.toVector());
  int direction = 6; // ==Down; this is right like this, because in TissueState::parseFromTrackedCellsFile a cell is started when a cell pixel was found *above* firstP
  do {
    if(p.isOnMargin()) {
      // cell has to be ignored 
      _ignoredCells.insert(id);
//      delete c; // delete cell again
//      return true;
    }
    std::vector<CellIndex> neighbors;
    direction = p.getLastCommonNeighboringBondPixelInCwOrientation(neighbors, id, direction);
    if(direction<0) {
      return false;
    }
    p.goTo(direction);
    // compute quantities
    Vector2D r(p.toVector());
    double areaTerm = lastR.x()*r.y() - lastR.y()*r.x();
    c->area += areaTerm;
    c->r += (r+lastR)*areaTerm;
    
    lastR = r;
  } while(p!=firstP);

  // finish computations
  c->area *= 0.5;
  c->r /= 6*c->area;

  // add cell to map
  _cells[id] = c;
  
  // **** second, compute cell elongation (needs computed cell center) ****
  // go around cell perimeter in ccw direction
  p = firstP;
  direction = 6; // ==Down; this is right like this, because in TissueState::parseFromTrackedCellsFile a cell is started when a cell pixel was found *above* firstP
  c->elongation.set(0,0);
  c->polarityR=c->polarityG=c->polarityB=c->elongation;
  c->intIntensityR=c->intIntensityG=c->intIntensityB=0;
  Vector2D r(p.toVector());
  Vector2D diff(r-c->r);
  double rLen=diff.norm(), rPhi=diff.angle();
  // for polarity:
  PixelValue oldData = OriginalFrame.data(p.position());
  unsigned char oldDataR=(oldData>>16)&0xFF, oldDataG=(oldData>>8)&0xFF, oldDataB=oldData&0xFF;
  do {
    std::vector<CellIndex> neighbors;
    direction = p.getLastCommonNeighboringBondPixelInCwOrientation(neighbors, id, direction);
    if(direction<0) {
      return false;
    }
    p.goTo(direction);

    // compute cell elongation
    double l = PixelFrame::NNeighborDistance[direction];
    Vector2D nextR(p.toVector());
    diff = nextR-c->r;
    double nextRLen=diff.norm(), nextRPhi=diff.angle();
    double deltaPhi = nextRPhi-rPhi;
    deltaPhi -= 2*M_PI*round(deltaPhi/2/M_PI);
    double rPhiMinusAlpha = M_PI + (nextR-r).angle();
    double sinAlpha = sin(rPhi-rPhiMinusAlpha);
    
    c->elongation += rLen*rLen*sinAlpha*sinAlpha*(
               nematicFromAngleAndNorm(rPhiMinusAlpha, l/(rLen*sinAlpha) - 2*deltaPhi)
              +nematicFromAngleAndNorm(rPhiMinusAlpha+0.25*M_PI, 2*log(rLen/nextRLen))
            );
    
    // polarity, separately for each color channel (R, G, B)
    PixelValue data = OriginalFrame.data(p.position());
    unsigned char dataR=(data>>16)&0xFF, dataG=(data>>8)&0xFF, dataB=data&0xFF;
    c->polarityR += nematicFromAngleAndNorm(rPhi+0.5*deltaPhi, deltaPhi*0.5*(dataR+oldDataR));
    c->intIntensityR += deltaPhi*0.5*(dataR+oldDataR);
    c->polarityG += nematicFromAngleAndNorm(rPhi+0.5*deltaPhi, deltaPhi*0.5*(dataG+oldDataG));
    c->intIntensityG += deltaPhi*0.5*(dataR+oldDataG);
    c->polarityB += nematicFromAngleAndNorm(rPhi+0.5*deltaPhi, deltaPhi*0.5*(dataB+oldDataB));
    c->intIntensityB += deltaPhi*0.5*(dataR+oldDataB);
    
    r = nextR;
    rLen = nextRLen;
    rPhi = nextRPhi;
    oldDataR = dataR;
    oldDataG = dataG;
    oldDataB = dataB;
  } while(p!=firstP);
  
  c->elongation /= 2*c->area;

  // **** third, go around again and create structures ****
  // go around cell perimeter
  int countFirstP = 0;
  p = firstP;
  std::vector<CellIndex> oldNeighbors;
  int oldDirection = -1;
  Vertex *curVertex = NULL;
  DirectedBond *curDirectedBond = NULL;
  double curBondLength = 0.0;
  do {
    // check neighboring cells and next direction
    std::vector<CellIndex> curNeighbors;
    int curDirection = p.getLastCommonNeighboringBondPixelInCwOrientation(curNeighbors, id, (oldDirection>=0)?oldDirection:6);
    if(curDirection<0) {
      return false;
    }
    
    bool isVertex = false;
    Pixel vertexP(p), enterP(p), leaveP(p);
    
    // -- check for thick vertex --
    std::set<Pixel> thickVertex;
    p.addToThickVertex(thickVertex);
    if(thickVertex.size()>1) {
//      if(thickVertex.size()>4) {
//        std::cout << "TissueState::addCell: thick vertex too big on pixel " << p.positionString() << ", direction " << oldDirection << ", for cell id " << std::hex << id << "." << std::endl;
//        return false;
//      }
      for(std::set<Pixel>::iterator it=thickVertex.begin(); it!=thickVertex.end(); ++it) {
        if(*it<vertexP) {
          vertexP = *it;
        }
      }

      isVertex = true;
      // go to last pixel end of thick vertex
      do {
        Pixel testP(p);
        testP.goTo(curDirection);
        if(thickVertex.count(testP)==0) {
          break;
        }
        p = testP;
//        curBondLength += PixelFrame::NNeighborDistance[curDirection]; // NOTE: thick vertices don't count in for bon length in order for both bonds of a conjugated pair to have the same length
        curDirection = p.getLastCommonNeighboringBondPixelInCwOrientation(curNeighbors, id, curDirection);
        if(curDirection<0) {
          return false;
        }
      } while(true);
      leaveP = p;
//      std::cout << "TissueState::addCell: thick vertex at " << p.positionString() << ", size " << thickVertex.size() << ", for cell id " << std::hex << id << std::dec << ", vertex " << vertexP.positionString() << ", old direction " << oldDirection << ", new direction " << curDirection <<  "." << std::endl;
    } else {
      // start with most probable case:
      if(curNeighbors.size()==2) {
        // consistency check
        if((oldNeighbors.size()>0) && (oldNeighbors.size()<3) && !compareCyclicLists(oldNeighbors, curNeighbors)) {
          std::cout << "TissueState::addCell: change of neighbor set, but no thick vertex at pixel " << p.positionString() << " for cell id " << std::hex << id << std::endl;
          return false;
        }
        // boundary between two cells
      } else if(curNeighbors.size()>2) {
        // vertex
        isVertex = true;
      } else if(curNeighbors.size()==1) {
        std::cout << "TissueState::addCell: Only one neighbor of pixel " << p.positionString() << " for cell id " << std::hex << id << std::endl;
        return false;
      } else if(curNeighbors.size()==0) {
        std::cout << "TissueState::addCell: No neighbors of pixel " << p.positionString() << " for cell id " << std::hex << id << std::endl;
        return false;
      }
    }

    if(isVertex) {
      // check, if come from valid direction
      if(oldDirection>=0) {
        // first, get pointer to vertex
        // check, if in map
        if(_vertexMap.count(vertexP)==0) {
          // no -> create
          curVertex = newVertex(vertexP.toVector());
          _vertexMap[vertexP] = curVertex;
        } else {
          curVertex = _vertexMap[vertexP];
        }

        // -- check out conjugated bond of previous bond --
        pair<Pixel, int> pixelAndDirectionConj(enterP, PixelFrame::oppositeDirection(oldDirection));
        DirectedBond *oldConjugatedBond = NULL;
        // old conjugated bond already there?
        if(_directedBondMap.count(pixelAndDirectionConj)>0) {
          oldConjugatedBond = _directedBondMap[pixelAndDirectionConj];
        } else {
          // no -> create
          oldConjugatedBond = newBond();
          _directedBondMap[pixelAndDirectionConj] = oldConjugatedBond;
        }

        // -- check out previous bond --
        // previous bond there?
        if(curDirectedBond) {
          // wire with cur vertex
          curDirectedBond->leftVertex = curVertex;
          // conjugated bond already there?
          if(oldConjugatedBond) {
            curDirectedBond->conjBond = oldConjugatedBond;
            oldConjugatedBond->conjBond = curDirectedBond;
          }
          // set length
          curDirectedBond->length = curBondLength;
        }

        // -- check out next bond --
        pair<Pixel, int> pixelAndDirection(leaveP, curDirection);
        if(_directedBondMap.count(pixelAndDirection)==0) {
          // new bond does not yet exist -> create
          curDirectedBond = newBond();
          _directedBondMap[pixelAndDirection] = curDirectedBond;
        } else {
          curDirectedBond = _directedBondMap[pixelAndDirection];
        }
        curBondLength = 0.0;
        
        // - add to cell -
        if(curDirectedBond->cell) {
          // cell already set?
          if(curDirectedBond->cell==c) {
            // right cell -> we're done!
            break;
          } else {
            // wrong cell
            std::cout << "TissueState::addCell: Bond already exists, but wrong cell at pixel " << p.positionString();
            std::cout << std::hex << " with cur cell id " << (*curNeighbors.begin()) << " and cell id in bond " << curDirectedBond->cell->id << std::endl;
            return false;
          }
        } else {
          // add
          curDirectedBond->cell = c;
          c->bonds.push_back(curDirectedBond);
        }

        // - add new bond (if appropriate) and (if appropriate) old conjugate bond to vertex -
        // have to add old conj. bond to vertex?
        if(!oldConjugatedBond->rightVertex) {
          // have to add new bond to vertex?
          if(!curDirectedBond->rightVertex) {
            // have to insert both into curVertex
            // first insert curBond at end
            curDirectedBond->rightVertex = curVertex;
            curVertex->bonds.insert(curVertex->bonds.end(), curDirectedBond);
            // then old conj bond
            oldConjugatedBond->rightVertex = curVertex;
            curVertex->bonds.insert(curVertex->bonds.end(), oldConjugatedBond);
          } else {
            // have to insert old conj. after new bond
            std::vector<DirectedBond *>::iterator it = std::find(curVertex->bonds.begin(), curVertex->bonds.end(), curDirectedBond);
            if(it!=curVertex->bonds.end()) {
              // found -> insert old conj. bond after -> in ccw direction
              oldConjugatedBond->rightVertex = curVertex;
              ++it;
              curVertex->bonds.insert(it, oldConjugatedBond);
            } else {
              // not found, although new bond was there, before
              std::cout << "TissueState::addCell: New bond was there before, but did not find it in vertex at pixel " << p.positionString() << " for cell id " << std::hex << c->id << std::endl;
              return false;
            }
          }
        } else {
          if(!curDirectedBond->rightVertex) {
            // find old conj. bond in vertex
            std::vector<DirectedBond *>::iterator it = std::find(curVertex->bonds.begin(), curVertex->bonds.end(), oldConjugatedBond);
            if(it!=curVertex->bonds.end()) {
              // found -> insert new bond before -> in cw direction
              curDirectedBond->rightVertex = curVertex;
              curVertex->bonds.insert(it, curDirectedBond);
            } else {
              // not found, although ocb was there, before
              std::cout << "TissueState::addCell: Old conjugated bond was there before, but did not find it in vertex at pixel " << p.positionString() << " for cell id " << std::hex << c->id << std::endl;
              return false;
            }
          } else {
            // both bonds are already connected!
          }
        }
      }
    }
    
    // next pixel
    p.goTo(curDirection);
    curBondLength += PixelFrame::NNeighborDistance[curDirection];
    if(p==firstP) {
      ++countFirstP;
      if(countFirstP>1) {
        // running in a circle, here
        // create artificial vertex and directed bond:
        curVertex = newVertex(firstP.toVector());
        curDirectedBond = newBond();
        curDirectedBond->rightVertex=curDirectedBond->leftVertex=curVertex;
        curDirectedBond->cell = c;
        c->bonds.push_back(curDirectedBond);
        curVertex->bonds.push_back(curDirectedBond);
        break;
      }
    }
    
    // keep old neighbors in mind
    oldNeighbors.swap(curNeighbors);
    // keep old directions in mind
    oldDirection = curDirection;
  } while(true);
  
  return true;
}

/** checks leftVertex<->conjBond consistency, bond sorting in vertices, and bond sorting in cells */
bool TissueState::checkTopologicalConsistency() const {
  const double Small = 1e-8;
  
  // leftVertex<->conjBond consistency and comparison of bond lengths
  for(TissueState::BondIterator it=beginBondIterator(); it!=endBondIterator(); ++it) {
    DirectedBond *b = bond(it);
    if(b->conjBond && (b->leftVertex!=b->conjBond->rightVertex)) {
      std::cout << "TissueState::checkTopologicalConsistency: inconsistency in bond between leftVertex = " << b->leftVertex << " and conjBond->rightVertex: " << b->conjBond->rightVertex << std::endl;
      std::cout << "TissueState::checkTopologicalConsistency: vertices " << b->rightVertex->r.x() << ", " << b->rightVertex->r.y() << " and " << b->conjBond->rightVertex->r.x() << ", " << b->conjBond->rightVertex->r.y() << std::endl;
      return false;
    }
    if(b->conjBond && (fabs(b->length-b->conjBond->length)>Small)) {
      std::cout << "TissueState::checkTopologicalConsistency: inconsistency in bond length " << b->rightVertex->r << "->" << b->leftVertex->r << "; diff: " << b->length-b->conjBond->length << std::endl;
    }
  }
  
  // vertex sorting
  for(unsigned int i=0; i<_vertices.size(); ++i) {
    Vertex *v = _vertices[i];
    if(v->bonds.size()>0) {
      Cell *lastCell = v->bonds[v->bonds.size()-1]->cell;
  //    std::cout << "TissueState::checkTopologicalConsistency: vertex " << v->r.x() << ", " << 1946-v->r.y() << ", size: " << v->bonds.size() << std::endl;
  //    std::cout << "TissueState::checkTopologicalConsistency: vertices " << v->bonds[v->bonds.size()-1]->rightVertex->r.x() << ", " << 1946-v->bonds[v->bonds.size()-1]->rightVertex->r.y() << " and " << v->bonds[v->bonds.size()-1]->leftVertex->r.x() << ", " << 1946-v->bonds[v->bonds.size()-1]->leftVertex->r.y() << std::endl;
      for(unsigned int j=0; j<v->bonds.size(); ++j) {
        DirectedBond *b = v->bonds[j];
  //      std::cout << "TissueState::checkTopologicalConsistency: vertices " << b->rightVertex->r.x() << ", " << 1946-b->rightVertex->r.y() << " and " << b->leftVertex->r.x() << ", " << 1946-b->leftVertex->r.y() << std::endl;
        if(b->conjBond && (lastCell!=b->conjBond->cell)) {
          std::cout << "TissueState::checkTopologicalConsistency: inconsistency in vertex sorting, vertex: " << v << std::endl;
          std::cout << "TissueState::checkTopologicalConsistency: inconsistency in vertex sorting, conjBond: " << b->conjBond << std::endl;
          std::cout << "TissueState::checkTopologicalConsistency: inconsistency in vertex sorting, conjBondCell: " << std::hex << b->conjBond->cell->id << std::endl;
          std::cout << "TissueState::checkTopologicalConsistency: inconsistency in vertex sorting, lastCell: " << lastCell->id << std::endl;
          std::cout << "TissueState::checkTopologicalConsistency: inconsistency in vertex sorting, bond: " << b << std::endl;
          std::cout << "TissueState::checkTopologicalConsistency: inconsistency in vertex sorting, cell: " << b->cell->id << std::endl;
  //        std::cout << "TissueState::checkTopologicalConsistency: cells " << std::hex << lastCell->id << ", " << b->conjBond->cell->id << "; " << b->cell->id << "." << std::endl;
  //         std::cout << "TissueState::checkTopologicalConsistency: vertices " << b->rightVertex->r.x() << ", " << 1946-b->rightVertex->r.y() << " and " << b->conjBond->rightVertex->r.x() << ", " << 1946-b->conjBond->rightVertex->r.y() << std::endl;
          return false;
        }
        lastCell = b->cell;
      }
    }
  }
  
  // cell sorting
  for(CellConstIterator it=beginCellIterator(); it!=endCellIterator(); ++it) {
    Cell *c = cell(it);
    if(c->bonds.size()>0) {
      Vertex *lastVertex = c->bonds[c->bonds.size()-1]->leftVertex;
      for(unsigned int j=0; j<c->bonds.size(); ++j) {
        DirectedBond *b = c->bonds[j];
        if(lastVertex!=b->rightVertex) {
          std::cout << "TissueState::checkTopologicalConsistency: inconsistency in cell sorting!" << std::endl;
          return false;
        }
        lastVertex = b->leftVertex;
      }
    }
  }
  
  return true;
}

void TissueState::markMarginVertices() {
  // reset
  for(unsigned int i=0; i<_vertices.size(); ++i) {
    _vertices[i]->margin = false;
  }
  // set
  for(TissueState::BondIterator it=beginBondIterator(); it!=endBondIterator(); ++it) {
    DirectedBond *b = bond(it);
    if(!b->conjBond) {
      b->rightVertex->margin = true;
    }
  }
  // consistency check
  for(TissueState::BondIterator it=beginBondIterator(); it!=endBondIterator(); ++it) {
    DirectedBond *b = bond(it);
    if(!b->conjBond && !b->leftVertex->margin) {
      std::cout << "TissueState::markMarginVertices: topological inconsistency!" << std::endl;
      throw std::exception();
    }
  }
}

void TissueState::removeMarginCells() {
  for(set<CellIndex>::iterator it=_ignoredCells.begin(); it!=_ignoredCells.end(); ++it) {
    if(contains(*it)) {
      removeCell(*it);
    } else {
      cout << "TissueState::removeMarginCells: did not find cell with id " << *it << "!" << endl;
      throw std::exception();
    }
  }
}


void recurseOverApoptoticNeighbors(Cell *c, std::set<Cell*> &connectedApoptoticCells, std::set<Cell*> &nonApoptoticNeighbors) {
  connectedApoptoticCells.insert(c);
  for(unsigned int j=0; j<c->bonds.size(); ++j) {
    if(!c->bonds[j]->conjBond) { // NOTE: can I lift this in principle?
      std::cout << "recurseOverApoptoticNeighbors: disappearing margin cell encountered!" << std::endl;
      throw std::exception();
    }
    Cell *n = c->bonds[j]->conjBond->cell;
    if((n->duringTransitionAfter!=Cell::Apoptosis) && (n->duringTransitionAfter!=Cell::SegmentationErrorDisappearance)) {
      // neighbor is not an apoptotic cell
      nonApoptoticNeighbors.insert(n);
    } else {
      // neighbor is an apoptotic cell
      if(!connectedApoptoticCells.count(n)) {
        recurseOverApoptoticNeighbors(n, connectedApoptoticCells, nonApoptoticNeighbors);
      }
    }
  }
}

TissueState *TissueState::createCopyFusingDivisions(const TissueState &previous) const {
  TissueState *s = new TissueState(_frameNumber, _time);
  
  // -- copy *this ---
  // maps
  std::map<Cell*,Cell*> oldCellToNewCell;
  std::map<Vertex*,Vertex*> oldVertexToNewVertex;
  std::map<DirectedBond*,DirectedBond*> oldBondToNewBond;
  // copy all cells
  for(CellConstIterator it=beginCellIterator(); it!=endCellIterator(); ++it) {
    Cell *c = cell(it);
    if(c->duringTransitionBefore==Cell::Divided) {
      Cell *sc = s->cellCheck(c->mother);
      if(!sc) {
        Cell *mother = previous.cellCheck(c->mother);
        if(!mother) {
          std::cout << "TissueState::createCopyFusingDivisions: mother cell not found in previous state!" << std::endl;
          throw std::exception();
        }
        Cell *sister = cellCheck(c->sister);
        sc = new Cell(*mother);
        s->_cells[mother->id] = sc;
        if(sister) {
          sc->r = (c->area*c->r + sister->area*sister->r)/(c->area + sister->area);
        } else {
          sc->r = c->r;
//          std::cout << "TissueState::createCopyFusingDivisions: sister cell not found in current state!" << std::endl;
//          throw std::exception();
        }
      }
      oldCellToNewCell[c] = sc;
    } else {
      Cell *pc = previous.cellCheck(c->id);
      Cell *sc;
      if(pc) {
        sc = new Cell(*pc);
      } else {
        sc = new Cell(c->id);
        sc->duringTransitionBefore = Cell::SegmentationErrorAppearance;
        sc->duringTransitionAfter = Cell::Stays;
      }
      sc->r = c->r;
      s->_cells[c->id] = sc;
      oldCellToNewCell[c] = sc;
    }
  }
  // copy all vertices
  for(int i=0; i<numberOfVertices(); ++i) {
    Vertex *v = vertex(i);
    Vertex *sv = new Vertex(*v);
    s->_vertices.push_back(sv);
    oldVertexToNewVertex[v] = sv;
  }
  // copy all bonds
  for(BondConstIterator it=beginBondIterator(); it!=endBondIterator(); ++it) {
    DirectedBond *b = bond(it);
    DirectedBond *sb = s->newBond();
    if(!b->cell) {
      std::cout << "TissueState::createCopyFusingDivisions: NULL cell pointer!" << std::endl;
      throw std::exception();
    }
    sb->cell = oldCellToNewCell[b->cell];
    if(!b->rightVertex) {
      std::cout << "TissueState::createCopyFusingDivisions: NULL right vertex pointer!" << std::endl;
      throw std::exception();
    }
    sb->rightVertex = oldVertexToNewVertex[b->rightVertex];
    if(!b->leftVertex) {
      std::cout << "TissueState::createCopyFusingDivisions: NULL left vertex pointer!" << std::endl;
      throw std::exception();
    }
    sb->leftVertex = oldVertexToNewVertex[b->leftVertex];
    oldBondToNewBond[b] = sb;
  }
  // wire bonds
  for(BondConstIterator it=beginBondIterator(); it!=endBondIterator(); ++it) {
    DirectedBond *b = bond(it);
    if(b->conjBond) {
      if(!(oldBondToNewBond[b]->conjBond = oldBondToNewBond[b->conjBond])) {
        std::cout << "TissueState::createCopyFusingDivisions: NULL new conj. bond pointer!" << std::endl;
        throw std::exception();
      }
    } else {
      oldBondToNewBond[b]->conjBond = NULL;
    }
  }
  // wire bonds within cells
  std::vector<DirectedBond*> newBondsToBeRemoved;
  std::map<Vertex*,pair<int,int> > oldDivisionVertices;
  std::set<Vertex *> newVerticesToBeChecked;
  for(CellConstIterator it=beginCellIterator(); it!=endCellIterator(); ++it) {
    Cell *c = cell(it);
    Cell *sc = oldCellToNewCell[c];
    if(c->duringTransitionBefore==Cell::Divided) {
      // find bond in c, the new version which connects to sc on both sides
      int jFound=-1;
      for(unsigned int j=0; j<c->bonds.size(); ++j) {
        DirectedBond *sbConj = oldBondToNewBond[c->bonds[j]]->conjBond;
        if(sbConj && (sbConj->cell==sc)) {
          jFound = j;
        }
      }
      if(jFound<0) {
        if(contains(c->sister)) {
//          std::cout << "TissueState::createCopyFusingDivisions: did not find division furrow bond within dividing cell, but sister cell is present!" << std::endl;
//          std::cout << _frameNumber << std::endl;
//          std::cout << std::hex << sc->id << " " << c->id << " " << cell(c->sister)->id << std::endl;
//          std::cout << std::dec << sc->r << " " << c->r << " " << cell(c->sister)->r << std::endl;
//          throw std::exception();
          // then mother and daughter might only be connected via a vertex; is this true?
          if(c->id==c->mother) {
            // look for vertex and indices
            Vertex *vFoundGlobal=NULL;
            int jFound=-1, kMFoundGlobal=-1, kDFoundGlobal=-1;
            for(unsigned int j=0; j<c->bonds.size(); ++j) {
              Vertex *v = c->bonds[j]->rightVertex;
              int kMFound = -1;
              for(unsigned int k=0; k<v->bonds.size(); ++k) {
                DirectedBond *b = v->bonds[k];
                if(b->cell->id==c->id) {
                  if(kMFound<0) {
                    kMFound = k;
                  } else {
                    std::cout << "TissueState::createCopyFusingDivisions: mother cell found twice around vertex!" << std::endl;
                    throw std::exception();
                  }
                }
                if(b->cell->id==c->sister) {
                  if(kDFoundGlobal<0) {
//                    std::cout << "TissueState::createCopyFusingDivisions: found: " << b << " " << b->cell << " " << v << std::endl;
                    kDFoundGlobal = k;
                  } else {
                    std::cout << "TissueState::createCopyFusingDivisions: daughter cell found twice around vertex or mother cell!" << std::endl;
                    throw std::exception();
                  }
                }
              }
              if((jFound<0) && (kDFoundGlobal>=0)) {
                jFound = j;
                kMFoundGlobal = kMFound;
                vFoundGlobal = v;
//                std::cout << "TissueState::createCopyFusingDivisions: vertex" << v << std::endl;
//                for(unsigned int k=0; k<v->bonds.size(); ++k) {
//                  DirectedBond *b = v->bonds[k];
////                  std::cout << "TissueState::createCopyFusingDivisions: bond " << b << "; cell " << std::hex << b->cell->id << "; conjBond " << b->conjBond << "; conjCell " << std::hex << b->conjBond->cell->id << std::endl;
//                  std::cout << "TissueState::createCopyFusingDivisions: new bond " << oldBondToNewBond[b] << "; new cell " << std::hex << oldBondToNewBond[b]->cell->id << "; new conjBond " << oldBondToNewBond[b]->conjBond << "; new conjCell " << std::hex << oldBondToNewBond[b]->conjBond->cell->id << std::endl;
//                }
              }
            }
            if((kDFoundGlobal<0) || (kMFoundGlobal<0)){
              std::cout << "TissueState::createCopyFusingDivisions: no vertex division found!" << std::endl;         
              throw std::exception();
            }
//            std::cout << "TissueState::createCopyFusingDivisions: vertex division: " << kMFoundGlobal << " " << kDFoundGlobal << " " << vFoundGlobal << std::endl;
            // store vertex and indices            
            oldDivisionVertices.insert(pair<Vertex*,pair<int,int> >(vFoundGlobal,pair<int,int>(kMFoundGlobal,kDFoundGlobal)));
            
            // add mother cell bonds in sc, start from division vertex
            for(unsigned int j=0; j<c->bonds.size(); ++j) {
              DirectedBond *sb = oldBondToNewBond[c->bonds[(j+jFound)%c->bonds.size()]];
              if(!sb) {
                std::cout << "TissueState::createCopyFusingDivisions: NULL new bond pointer within dividing cell!" << std::endl;
                throw std::exception();
              }
              sc->bonds.push_back(sb);
            }
            // look for vertex in daugther cell
            Cell *d = cell(c->sister);
//            std::cout << "TissueState::createCopyFusingDivisions: " << std::hex << d->id << " " << d << std::endl;
            jFound=-1;
            for(unsigned int j=0; j<d->bonds.size(); ++j) {
//              std::cout << "TissueState::createCopyFusingDivisions: " << d->bonds[j] << " " << d->bonds[j]->rightVertex << std::endl;
              if(d->bonds[j]->rightVertex==vFoundGlobal) {
                jFound = j;
              }
            }
            if(jFound<0) {
              std::cout << "TissueState::createCopyFusingDivisions: vertex not found around daughter cell!" << std::endl;
              throw std::exception();
            }
            // add daughter cell bonds in sc, start from division vertex
            for(unsigned int j=0; j<d->bonds.size(); ++j) {
              DirectedBond *sb = oldBondToNewBond[d->bonds[(j+jFound)%d->bonds.size()]];
              if(!sb) {
                std::cout << "TissueState::createCopyFusingDivisions: NULL new bond pointer within dividing cell!" << std::endl;
                throw std::exception();
              }
              sc->bonds.push_back(sb);
            }
          }
        } else {
          // then, sister cell of c may not be part of the tissue state, here
          // just copy bonds
          for(unsigned int j=0; j<c->bonds.size(); ++j) {
            DirectedBond *sb = oldBondToNewBond[c->bonds[j]];
            if(!sb) {
              std::cout << "TissueState::createCopyFusingDivisions: NULL new bond pointer within cell!" << std::endl;
              throw std::exception();
            }
            sc->bonds.push_back(sb);
          }
        }
      } else {
        // add new bond to bonds to be removed
        newBondsToBeRemoved.push_back(oldBondToNewBond[c->bonds[jFound]]);
        // create neighbors in sc, start from division furrow but leave it out
        for(unsigned int j=1; j<c->bonds.size(); ++j) {
          DirectedBond *sb = oldBondToNewBond[c->bonds[(j+jFound)%c->bonds.size()]];
          if(!sb) {
            std::cout << "TissueState::createCopyFusingDivisions: NULL new bond pointer within dividing cell!" << std::endl;
            throw std::exception();
          }
          sc->bonds.push_back(sb);
        }
      }
    } else {
      // just copy bonds
      for(unsigned int j=0; j<c->bonds.size(); ++j) {
        DirectedBond *sb = oldBondToNewBond[c->bonds[j]];
        if(!sb) {
          std::cout << "TissueState::createCopyFusingDivisions: NULL new bond pointer within cell!" << std::endl;
          throw std::exception();
        }
        sc->bonds.push_back(sb);
      }
    }
  }
  // wire bonds within vertices
  for(int i=0; i<numberOfVertices(); ++i) {
    Vertex *v = vertex(i);
    // is v a division vertex?
    if(oldDivisionVertices.count(v)) {
      // --- spit vertex ---
      // get split indices
      int kM = oldDivisionVertices[v].first, kD = oldDivisionVertices[v].second;
      
      // already created vertex
      Vertex *sv = oldVertexToNewVertex[v];
      
      // create second vertex:
      Vertex *svSecond = new Vertex(*v);
      s->_vertices.push_back(svSecond);
      
      // add bonds
      int count = 0;
      for(unsigned int k=0; k<v->bonds.size(); ++k) {
        DirectedBond *sb = oldBondToNewBond[v->bonds[k]];
        if(!sb) {
          std::cout << "TissueState::createCopyFusingDivisions: NULL deformed bond pointer within vertex!" << std::endl;
          throw std::exception();
        }
        if(count%2==0) {
          sv->bonds.push_back(sb);
        } else {
          svSecond->bonds.push_back(sb);
          sb->rightVertex = svSecond;
          if(sb->conjBond) {
            sb->conjBond->leftVertex = svSecond;
          }
        }
        if((k==kM) || (k==kD)) {
          ++count;
        }
      }
      
//      std::cout << "TissueState::createCopyFusingDivisions: splitted vertex " << v << " into " << sv << " and " << svSecond << "." << std::endl;
//      std::cout << "TissueState::createCopyFusingDivisions: neighbors: " << v->bonds.size() << ", " << sv ->bonds.size() << ", and " << svSecond->bonds.size() << "." << std::endl;
      newVerticesToBeChecked.insert(sv);
      newVerticesToBeChecked.insert(svSecond);
    } else { // no division vertex
      Vertex *sv = oldVertexToNewVertex[v];
      for(unsigned int j=0; j<v->bonds.size(); ++j) {
        DirectedBond *sb = oldBondToNewBond[v->bonds[j]];
        if(!sb) {
          std::cout << "TissueState::createCopyFusingDivisions: NULL deformed bond pointer within vertex!" << std::endl;
          throw std::exception();
        }
        sv->bonds.push_back(sb);
      }
    }
  }
  // remove bonds
  for(unsigned int i=0; i<newBondsToBeRemoved.size(); ++i) {
    newVerticesToBeChecked.insert(newBondsToBeRemoved[i]->rightVertex);
    s->removeBondWithoutCell(newBondsToBeRemoved[i]);
  }
  // check vertices
  for(std::set<Vertex*>::iterator it=newVerticesToBeChecked.begin(); it!=newVerticesToBeChecked.end(); ++it) {
    Vertex *v = *it;
    if(v->bonds.size()<=2) {
      // count actual number of bonds
      int actualNumberOfBonds = v->bonds.size();
      for(unsigned int j=0; j<v->bonds.size(); ++j) {
        if(!v->bonds[j]->conjBond) {
          ++actualNumberOfBonds;
        }
      }
      // check
      if(actualNumberOfBonds==2) {
        if(v->bonds.size()==2) {
          DirectedBond *b0 = v->bonds[0], *b1 = v->bonds[1];
          // check for special cases
          if((b0->conjBond==b1) || (b1->conjBond==b0)) {
            std::cout << "TissueState::createCopyFusingDivisions: strange 2 fold vertex in frame " << _frameNumber << " at " << v->r << "!" << std::endl;
            throw std::exception();
          }
          // connect 
          b0->conjBond->conjBond = b1->conjBond;
          b1->conjBond->conjBond = b0->conjBond;
          b0->conjBond->leftVertex = b1->leftVertex;
          b1->conjBond->leftVertex = b0->leftVertex;
          // remove bonds -> vertex is also removed
          b0->conjBond = NULL;
          b1->conjBond = NULL;
          s->removeBond(b0);
          s->removeBond(b1);
        } else if(v->bonds.size()==1) {
          DirectedBond *b0 = v->bonds[0];
          // look for other bond
          DirectedBond *b1 = NULL;
          for(BondIterator it=s->beginBondIterator(); it!=s->endBondIterator(); ++it) {
            if(bond(it)->leftVertex==v) {
              b1 = bond(it);
            }
          }
          if(!b1) {
            std::cout << "TissueState::createCopyFusingDivisions: margin vertex: other bond not found!" << std::endl;
            throw std::exception();
          }
          if(b0->conjBond || b1->conjBond) {
            std::cout << "TissueState::createCopyFusingDivisions: margin vertex: inconsistency!" << std::endl;
            throw std::exception();
          }
          // connect 
          b1->leftVertex = b0->leftVertex;
          // remove bond -> vertex is also removed
          s->removeBond(b0);
        } else {
          std::cout << "TissueState::createCopyFusingDivisions: Vertex has 2 actual bonds, but no direct bonds!" << std::endl;
          throw std::exception();
        }
      } 
      else if(actualNumberOfBonds<2) {
        std::cout << "TissueState::createCopyFusingDivisions: " << v->bonds.size() << " fold vertex!" << std::endl;
        throw std::exception();
      }
    }
  } 

  return s;
}


TissueState *TissueState::createCopyAndMoveCellPositionsTo(const TissueState &other) {
  TissueState *s = new TissueState(_frameNumber, _time);
  
  // -- copy *this ---
  // copy all cells
  for(CellConstIterator it=beginCellIterator(); it!=endCellIterator(); ++it) {
    Cell *c = cell(it);
    Cell *sc = new Cell(*c);
    s->_cells[c->id] = sc;
    c->deformedCell = sc;
  }
  // copy all vertices
  for(int i=0; i<numberOfVertices(); ++i) {
    Vertex *v = vertex(i);
    Vertex *sv = new Vertex(*v);
    s->_vertices.push_back(sv);
    v->deformedVertex = sv;
  }
  // copy all bonds
  for(BondConstIterator it=beginBondIterator(); it!=endBondIterator(); ++it) {
    DirectedBond *b = bond(it);
    DirectedBond *sb = s->newBond();
    b->deformedBond = sb;
    if(!b->cell || !(sb->cell = b->cell->deformedCell)) {
      std::cout << "TissueState::createCopyAndMoveCellPositionsTo: NULL cell pointer!" << std::endl;
      throw std::exception();
    }
    if(!b->rightVertex || !(sb->rightVertex = b->rightVertex->deformedVertex)) {
      std::cout << "TissueState::createCopyAndMoveCellPositionsTo: NULL right vertex pointer!" << std::endl;
      throw std::exception();
    }
    if(!b->leftVertex || !(sb->leftVertex = b->leftVertex->deformedVertex)) {
      std::cout << "TissueState::createCopyAndMoveCellPositionsTo: NULL left vertex pointer!" << std::endl;
      throw std::exception();
    }
  }
  // wire bonds
  for(BondConstIterator it=beginBondIterator(); it!=endBondIterator(); ++it) {
    DirectedBond *b = bond(it);
    if(b->conjBond) {
      if(!(b->deformedBond->conjBond = b->conjBond->deformedBond)) {
        std::cout << "TissueState::createCopyAndMoveCellPositionsTo: NULL deformed conj. bond pointer!" << std::endl;
        throw std::exception();
      }
    } else {
      b->deformedBond->conjBond = NULL;
    }
  }
  // wire bonds within cells
  for(CellConstIterator it=beginCellIterator(); it!=endCellIterator(); ++it) {
    Cell *c = cell(it);
    Cell *sc = c->deformedCell;
    for(unsigned int j=0; j<c->bonds.size(); ++j) {
      if(!c->bonds[j]->deformedBond) {
        std::cout << "TissueState::createCopyAndMoveCellPositionsTo: NULL deformed bond pointer within cell!" << std::endl;
        throw std::exception();
      }
      sc->bonds.push_back(c->bonds[j]->deformedBond);
    }
  }
  // wire bonds within vertices
  for(int i=0; i<numberOfVertices(); ++i) {
    Vertex *v = vertex(i);
    Vertex *sv = v->deformedVertex;
    for(unsigned int j=0; j<v->bonds.size(); ++j) {
      if(!v->bonds[j]->deformedBond) {
        std::cout << "TissueState::createCopyAndMoveCellPositionsTo: NULL deformed bond pointer within vertex!" << std::endl;
        throw std::exception();
      }
      sv->bonds.push_back(v->bonds[j]->deformedBond);
    }
  }

  // --- copy cell positions ---
  std::set<Cell*> affectedApoptoticCells;
  for(CellConstIterator it=s->beginCellIterator(); it!=s->endCellIterator(); ++it) {
    Cell *sc = s->cell(it);
    Cell *oc = other.cellCheck(sc->id);
    if(oc) {
      // set to position in the other state
      sc->r = oc->r;
    } else {
      // avoid unnecessary repetitions
      if(affectedApoptoticCells.count(sc)) {
        continue;
      }
      affectedApoptoticCells.insert(sc);
      // cell not there in other state
      if((sc->duringTransitionAfter!=Cell::Apoptosis) && (sc->duringTransitionAfter!=Cell::SegmentationErrorDisappearance)) {
        std::cout << "TissueState::createCopyAndMoveCellPositionsTo: frame " << _frameNumber << ": cell " << sc->id << sc->r << " not found in \'other\' state, but no apoptosis or segmentation error!" << std::endl;
        throw std::exception();
      }
      // get all connected apoptotic cells
      std::set<Cell*> connectedApoptoticCells, nonApoptoticNeighbors;
      recurseOverApoptoticNeighbors(sc, connectedApoptoticCells, nonApoptoticNeighbors);
      if(nonApoptoticNeighbors.size()==0) {
        std::cout << "TissueState::createCopyAndMoveCellPositionsTo: no non-apoptotic neighbor found!" << std::endl;
        throw std::exception();
      }
      // get the center of the surrounding cells in the other state
      Vector2D center(0,0);
      for(std::set<Cell*>::iterator it=nonApoptoticNeighbors.begin(); it!=nonApoptoticNeighbors.end(); ++it) {
        Cell *oc = other.cellCheck((*it)->id);
        if(!oc) {
          std::cout << "TissueState::createCopyAndMoveCellPositionsTo: cell not found in \'other\' state, but no apoptosis or segmentation error!" << std::endl;
          throw std::exception();
        }
        center += oc->r;
      }
      center /= nonApoptoticNeighbors.size();
      // set the positions of the apoptotic cells
      for(std::set<Cell*>::iterator it=connectedApoptoticCells.begin(); it!=connectedApoptoticCells.end(); ++it) {
        (*it)->r = center;
      }
    }
  }
  
  _deformedState = s;
  return s;
}

void TissueState::createTriangles() {
  for(int i=0; i<numberOfVertices(); ++i) {
    vertex(i)->createTriangles();
  }
  for(TissueState::BondIterator it=beginBondIterator(); it!=endBondIterator(); ++it) {
    bond(it)->createTriangle();
  }
}

double TissueState::totalAreaOfDualMargin() const {
  double twoTimesTotalArea = 0.0;
  for(int i=0; i<numberOfVertices(); ++i) {
    Vertex *v = vertex(i);
    if(v->margin) {
      for(unsigned int j=0; j<v->bonds.size(); ++j) {
        DirectedBond *b = v->bonds[j];
        if(b->conjBond) {
          Vector2D &r0 = b->cell->r, &r1 = b->conjBond->cell->r;
          // r1 lies in ccw direction around the margin from r0
          twoTimesTotalArea += -r1.x()*r0.y() + r0.x()*r1.y();
        }
      }
    }
  }
  return 0.5*twoTimesTotalArea;  
}

Matrix2x2 TissueState::totalDeformationOfDualMargin() const {
  Matrix2x2 sumMatrix(0,0,0,0);
  for(int i=0; i<numberOfVertices(); ++i) {
    Vertex *v = vertex(i);
    if(v->margin) {
      for(unsigned int j=0; j<v->bonds.size(); ++j) {
        DirectedBond *b = v->bonds[j];
        if(b->conjBond) {
          Vector2D &r0 = b->cell->r, &r1 = b->conjBond->cell->r;
          Vector2D &deformedR0 = b->deformedBond->cell->r, &deformedR1 = b->deformedBond->conjBond->cell->r;
          Vector2D deltaR(r1-r0);
          Vector2D rotDeltaR(deltaR.y(), -deltaR.x()); // rotated by pi/2 in cw sense -> normal vector pointing outwards, weighted with length
          sumMatrix += dyadicProduct(rotDeltaR, (deformedR1 - r1) + (deformedR0 - r0));
        }
      }
    }
  }  
  return sumMatrix/(2*totalAreaOfDualMargin());  
}

void TissueState::drawTime(QtImage &img, const Color &col, const double TextHeight) const {
  const double Factor=5.0; //, TextHeight=img.height()*heightPortion;
  char timeStr[1024];
  sprintf(timeStr, "%.1fhAPF", time());
  img.drawRectD(img.width()-(Factor+0.5)*TextHeight, 0, (Factor+0.5)*TextHeight, 1.7*TextHeight, Black);
  img.drawText(img.width()-Factor*TextHeight, 0.5*TextHeight, TextHeight, col, timeStr);
}

void TissueState::drawMargin(QtImage &img, const IsotropicTransformation &Ref2Pixel, const Color &col, const double thickness) const {
  for(TissueState::BondIterator it=beginBondIterator(); it!=endBondIterator(); ++it) {
    DirectedBond *b = bond(it);
    if(!b->conjBond) {
      img.drawBar(Ref2Pixel.map(b->rightVertex->r), Ref2Pixel.map(b->leftVertex->r), col, thickness);
    }
  }
}
  
void TissueState::drawBonds(QtImage &img, const IsotropicTransformation &Ref2Pixel, const Color &col, const double thickness, const bool includeMargin) const {
  std::set<DirectedBond*> drawnBonds;
  for(TissueState::BondIterator it=beginBondIterator(); it!=endBondIterator(); ++it) {
    DirectedBond *b = bond(it);
    if((includeMargin || b->conjBond) && !drawnBonds.count(b->conjBond)) {
      img.drawBar(Ref2Pixel.map(b->rightVertex->r), Ref2Pixel.map(b->leftVertex->r), col, thickness);
      drawnBonds.insert(b);
    }
  }
}

void TissueState::drawVertices(QtImage &img, const IsotropicTransformation &Ref2Pixel, const Color &col, const double radius) const {
  for(int i=0; i<numberOfVertices(); ++i) {
    Vertex *v = vertex(i);
    //DEBUG
    if(v->bonds.size()<2) {
      img.drawCircle(Ref2Pixel.map(v->r), radius, Blue);
    } else {
      img.drawCircle(Ref2Pixel.map(v->r), radius, col);
    }
  }  
}

void TissueState::drawCells(QtImage &img, const IsotropicTransformation &Ref2Pixel, const Color &col, const double radius) const {
  for(CellConstIterator it=beginCellIterator(); it!=endCellIterator(); ++it) {
    Cell *c = cell(it);
    img.drawCircle(Ref2Pixel.map(c->r), radius, col);
  }
}

void TissueState::drawDivisionPairs(QtImage &img, const IsotropicTransformation &Ref2Pixel, const Color &motherCol, const Color &daughterCol, const double radius, const Color &connectionCol, const double thickness) const {
  for(CellConstIterator it=beginCellIterator(); it!=endCellIterator(); ++it) {
    Cell *c = cell(it);
    if((c->duringTransitionBefore==Cell::Divided) && (c->mother==c->id)) {
      img.drawCircle(Ref2Pixel.map(c->r), radius, motherCol);
      Cell *oc = cell(c->sister);
      img.drawBar(Ref2Pixel.map(c->r), Ref2Pixel.map(oc->r), connectionCol, thickness);
    }
    if((c->duringTransitionBefore==Cell::Divided) && (c->mother!=c->id)) {
      img.drawCircle(Ref2Pixel.map(c->r), radius, daughterCol);
    }
  }
}

void TissueState::drawCellsTransitionBefore(QtImage &img, const IsotropicTransformation &Ref2Pixel, const double radius) const {
  for(CellConstIterator it=beginCellIterator(); it!=endCellIterator(); ++it) {
    Cell *c = cell(it);
    Color col(Blue);
    switch(c->duringTransitionBefore) {
      case Cell::UnclassifiedBefore:
        col = White;
        break;
      case Cell::Stayed:
        col = Blue;
        break;
      case Cell::Divided:
        col = Green;
        break;
      case Cell::MovedIntoMask:
        col = Yellow;
        break;
      case Cell::SegmentationErrorAppearance:
        col = Magenta;
        break;
      default:
        col = Black;
        break;
    }
    img.drawCircle(Ref2Pixel.map(c->r), radius, col);
  }
}

void TissueState::drawCellsTransitionAfter(QtImage &img, const IsotropicTransformation &Ref2Pixel, const double radius) const {
  for(CellConstIterator it=beginCellIterator(); it!=endCellIterator(); ++it) {
    Cell *c = cell(it);
    Color col(Blue);
    switch(c->duringTransitionAfter) {
      case Cell::UnclassifiedAfter:
        col = White;
        break;
      case Cell::Stays:
        col = Blue;
        break;
      case Cell::Divides:
        col = Green;
        break;
      case Cell::Apoptosis:
        col = Red;
        break;
      case Cell::MovesOutOfMask:
        col = Yellow;
        break;
      case Cell::SegmentationErrorDisappearance:
        col = Magenta;
        break;
      default:
        col = Black;
        break;
    }
    img.drawCircle(Ref2Pixel.map(c->r), radius, col);
  }
}

void TissueState::drawCellsPolarityNormalizedByIntIntensityR(QtImage &img, const IsotropicTransformation &Ref2Pixel, const Color &col, const double maxLenInPixels, const double width) const {
  for(CellConstIterator it=beginCellIterator(); it!=endCellIterator(); ++it) {
    Cell *c = cell(it);
    img.drawBar(Ref2Pixel.map(c->r), rotatedObject(maxLenInPixels*c->polarityR/c->intIntensityR, Ref2Pixel.angle()), col, width);
  }
}

#ifdef USE_NETCDF

bool TissueState::load(const std::string Filename) {
  cleanUp();
  
  // Change the error behavior of the netCDF C++ API
  NcError err(NcError::verbose_nonfatal);
  // Open the file.
  NcFile dataFile(Filename.c_str(), NcFile::ReadOnly);
  // Check to see if the file was opened
  if(!dataFile.is_valid()) return false;

  // load tissue state data
  NcVar *frameNumberVar, *timeVar;
  if(!(frameNumberVar = dataFile.get_var("TissueState::frameNumber"))) return false;
  _frameNumber = frameNumberVar->as_int(0);
  if(!(timeVar = dataFile.get_var("TissueState::time"))) return false;
  _time = timeVar->as_double(0);

  // load vertices
  Vertex::DataFileIO vio;
  if(!vio.load(*this, dataFile)) {
    return false;
  }
  // load cells
  Cell::DataFileIO cio;
  if(!cio.load(*this, dataFile)) {
    return false;
  }
  // load bonds
  DirectedBond::DataFileIO bio;
  if(!bio.load(*this, dataFile)) {
    return false;
  }
  
  // check for topological inconsistencies
  if(!checkTopologicalConsistency()) {
    std::cout << "Topological inconsistencies in " << Filename << "!" << std::endl;
    cleanUp();
    return false;
  }
  markMarginVertices();
  
  return true;
}


bool TissueState::save(const std::string Filename) const {
  // Change the error behavior of the netCDF C++ API
  NcError err(NcError::verbose_fatal);
  // Create the file.
  NcFile dataFile(Filename.c_str(), NcFile::Replace);
  // Check to see if the file was created.
  if(!dataFile.is_valid()) return false;
  
  // save tissue state data
  NcVar *frameNumberVar, *timeVar;
  if(!(frameNumberVar = dataFile.add_var("TissueState::frameNumber", ncInt))) return false;
  if(!frameNumberVar->put(&_frameNumber)) return false;
  if(!(timeVar = dataFile.add_var("TissueState::time", ncDouble))) return false;
  if(!timeVar->put(&_time)) return false;
  
  // save vertices
  Vertex::DataFileIO vio;
  if(!vio.save(dataFile, *this)) {
    return false;
  }
  // save cells
  Cell::DataFileIO cio;
  if(!cio.save(dataFile, *this)) {
    return false;
  }
  // save bonds
  DirectedBond::DataFileIO bio;
  if(!bio.save(dataFile, *this)) {
    return false;
  }
  
  return true;
}

#endif  

bool TissueState::exportToDbTables(const int ImageHeight, ofstream& framesFile, ofstream &verticesFile, ofstream &cellsFile, ofstream &ignoredCellsFile, ofstream &undirectedBondsFile, ofstream &directedBondsFile, unsigned long &lastVid, unsigned long &lastDbid, unsigned long &lastUbid) const {
  if(contains(Cell::VoidCellId)) {
    std::cout << "Found void cell id " << Cell::VoidCellId << " in frame " << _frameNumber << "!" << std::endl;
    return false;
  }
  
  // add data to frames table
  framesFile << _frameNumber << LogFile::DataSeparator << _time << '\n';
  if(framesFile.bad()) {
    std::cout << "Error writing frame data!" << std::endl;
    return false;
  }
  
  // add data to vertex table and generate map: vertex pointer -> vid
  std::map<Vertex*,unsigned long> vertexPointerToVid;
  for(int i=0; i<numberOfVertices(); ++i) {
    Vertex *v = vertex(i);
    ++lastVid;
    vertexPointerToVid[v] = lastVid;
    // write data:
    verticesFile << _frameNumber << LogFile::DataSeparator << lastVid << LogFile::DataSeparator << v->r.x() << LogFile::DataSeparator << ImageHeight-1-v->r.y() << '\n';
    if(verticesFile.bad()) {
      std::cout << "Error writing vertex data!" << std::endl;
      return false;
    }
  }
  
  // add data to cell table and generate map for bond sorting within cells
  std::map<DirectedBond*,DirectedBond*> leftBondSeenFromCellMap;
  // void cell
  cellsFile << _frameNumber << LogFile::DataSeparator << Cell::VoidCellId
          << LogFile::DataSeparator << 0 << LogFile::DataSeparator << 0 << LogFile::DataSeparator << 0 
          << LogFile::DataSeparator << -1 << LogFile::DataSeparator << -1 << LogFile::DataSeparator << 0 
          << LogFile::DataSeparator << 0 << LogFile::DataSeparator << 0 
          << LogFile::DataSeparator << 0 << LogFile::DataSeparator << 0 << LogFile::DataSeparator << 0 
          << LogFile::DataSeparator << 0 << LogFile::DataSeparator << 0 << LogFile::DataSeparator << 0 
          << LogFile::DataSeparator << 0 << LogFile::DataSeparator << 0 << LogFile::DataSeparator << 0 
          << '\n';
  // other cells
  for(TissueState::CellConstIterator it=beginCellIterator(); it!=endCellIterator(); ++it) {
    Cell *c = cell(it);
    for(unsigned int j=0; j<c->bonds.size(); ++j) {
      leftBondSeenFromCellMap[c->bonds[j]] = c->bonds[(j+1)%c->bonds.size()];
    }
    // write data:
    cellsFile << _frameNumber << LogFile::DataSeparator << c->id
            << LogFile::DataSeparator << (int)c->duringTransitionBefore << LogFile::DataSeparator << (int)c->duringTransitionAfter << LogFile::DataSeparator << c->daughter
            << LogFile::DataSeparator << c->r.x() << LogFile::DataSeparator << ImageHeight-1-c->r.y() << LogFile::DataSeparator << c->area
            << LogFile::DataSeparator << c->elongation.c1() << LogFile::DataSeparator << c->elongation.c2()
            << LogFile::DataSeparator << c->polarityR.c1() << LogFile::DataSeparator << c->polarityR.c2() << LogFile::DataSeparator << c->intIntensityR
            << LogFile::DataSeparator << c->polarityG.c1() << LogFile::DataSeparator << c->polarityG.c2() << LogFile::DataSeparator << c->intIntensityG
            << LogFile::DataSeparator << c->polarityB.c1() << LogFile::DataSeparator << c->polarityB.c2() << LogFile::DataSeparator << c->intIntensityB
            << '\n';
    if(cellsFile.bad()) {
      std::cout << "Error writing cell data!" << std::endl;
      return false;
    }
  }

  // ignored cells
  for(set<CellIndex>::const_iterator it=_ignoredCells.begin(); it!=_ignoredCells.end(); ++it) {
    // write data:
    ignoredCellsFile << _frameNumber << LogFile::DataSeparator << *it << '\n';
    if(ignoredCellsFile.bad()) {
      std::cout << "Error writing data for ignored cells!" << std::endl;
      return false;
    }
  }

   
  // add data to undirected bonds table and generate map: bond pointer -> dbid
  std::map<DirectedBond*,unsigned long> directedBondPointerToDbid;
  std::map<DirectedBond*,unsigned long> directedBondPointerToUbid;
  std::map<DirectedBond*,unsigned long> newConjBondDbids;
  std::map<DirectedBond*,DirectedBond*> newConjBondLeftOfConjBondAsSeenFromVoidCell;
  for(TissueState::BondConstIterator it=beginBondIterator(); it!=endBondIterator(); ++it) {
    DirectedBond *b = bond(it);
    ++lastDbid;
    directedBondPointerToDbid[b] = lastDbid;
    if(!b->conjBond) {
      ++lastDbid;
      newConjBondDbids[b] = lastDbid;
      
      // conj bond of bond left of conj bond seen from void cell
      Vertex *v = b->rightVertex;
      unsigned int i;
      for(i=0; i<v->bonds.size(); ++i) {
        if(v->bonds[i]==b) break;
      }
      if(i==v->bonds.size()) {
        std::cout << "TissueState::exportToDbTables: bond not found within rightVertex!" << std::endl;
        throw std::exception();
      }
      DirectedBond *oneBondCwAtVertex = v->bonds[ (i+v->bonds.size()-1)%v->bonds.size() ];
      
      Cell *c = oneBondCwAtVertex->cell;
      for(i=0; i<c->bonds.size(); ++i) {
        if(c->bonds[i]==oneBondCwAtVertex) break;
      }
      if(i==c->bonds.size()) {
        std::cout << "TissueState::exportToDbTables: bond not found within cell!" << std::endl;
        throw std::exception();
      }
      DirectedBond *oneBondCwAtCell = c->bonds[ (i+c->bonds.size()-1)%c->bonds.size() ];
      
      newConjBondLeftOfConjBondAsSeenFromVoidCell[b] = oneBondCwAtCell;
    }
    if((!b->conjBond) || (directedBondPointerToUbid.count(b->conjBond)==0)) {
      ++lastUbid;
      directedBondPointerToUbid[b] = lastUbid;
      undirectedBondsFile << _frameNumber << LogFile::DataSeparator << lastUbid << LogFile::DataSeparator << b->length << '\n';
      if(undirectedBondsFile.bad()) {
        std::cout << "Error writing undirectedBond data!" << std::endl;
        return false;
      }     
    } else {
      directedBondPointerToUbid[b] = directedBondPointerToUbid.at(b->conjBond);
    }
  }
  
  // add data to directed bonds table
  for(TissueState::BondConstIterator it=beginBondIterator(); it!=endBondIterator(); ++it) {
    DirectedBond *b = bond(it);
    if(b->conjBond) {
      directedBondsFile << _frameNumber << LogFile::DataSeparator << directedBondPointerToDbid.at(b) << LogFile::DataSeparator << directedBondPointerToDbid.at(b->conjBond) << LogFile::DataSeparator << directedBondPointerToUbid.at(b) << LogFile::DataSeparator << b->cell->id << LogFile::DataSeparator << vertexPointerToVid[b->rightVertex] << LogFile::DataSeparator << directedBondPointerToDbid.at(leftBondSeenFromCellMap.at(b)) << '\n';
      if(directedBondsFile.bad()) {
        std::cout << "Error writing directedBond data!" << std::endl;
        return false;
      }     
    } else {
      // add bond
      directedBondsFile << _frameNumber << LogFile::DataSeparator << directedBondPointerToDbid.at(b) << LogFile::DataSeparator << newConjBondDbids.at(b) << LogFile::DataSeparator << directedBondPointerToUbid.at(b) << LogFile::DataSeparator << b->cell->id << LogFile::DataSeparator << vertexPointerToVid[b->rightVertex] << LogFile::DataSeparator << directedBondPointerToDbid.at(leftBondSeenFromCellMap.at(b)) << '\n';
      if(directedBondsFile.bad()) {
        std::cout << "Error writing directedBond data!" << std::endl;
        return false;
      }
      // add new conj bond
      DirectedBond *conjBondOfLeftBondAsSeenFromVoidCell = newConjBondLeftOfConjBondAsSeenFromVoidCell.at(b);
      if(newConjBondDbids.count(conjBondOfLeftBondAsSeenFromVoidCell)==0) {
        std::cout << "TissueState::exportToDbTables: Problem with margin bonds!" << std::endl;
        throw std::exception();
      }
      directedBondsFile << _frameNumber << LogFile::DataSeparator << newConjBondDbids.at(b) << LogFile::DataSeparator << directedBondPointerToDbid.at(b) << LogFile::DataSeparator << directedBondPointerToUbid.at(b) << LogFile::DataSeparator << Cell::VoidCellId << LogFile::DataSeparator << vertexPointerToVid[b->leftVertex] << LogFile::DataSeparator << newConjBondDbids.at(conjBondOfLeftBondAsSeenFromVoidCell) << '\n';
      if(directedBondsFile.bad()) {
        std::cout << "Error writing directedBond data!" << std::endl;
        return false;
      }     
    }
  }
  
  return true;
}
