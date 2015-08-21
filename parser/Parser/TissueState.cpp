#include <set>
#include <stack>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>

#include "TextProgressBar.h"
#include "QtRasterImage.h"

#include "LogFile.h"
#include "Pixel.h"
#include "TissueState.h"

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
    // write data; revert signs of nematic xy components such that the values in the .dat files corresponds to a coordinate system with the y axis pointing upwards
    cellsFile << _frameNumber << LogFile::DataSeparator << c->id
            << LogFile::DataSeparator << (int)c->duringTransitionBefore << LogFile::DataSeparator << (int)c->duringTransitionAfter << LogFile::DataSeparator << c->daughter
            << LogFile::DataSeparator << c->r.x() << LogFile::DataSeparator << ImageHeight-1-c->r.y() << LogFile::DataSeparator << c->area
            << LogFile::DataSeparator << c->elongation.c1() << LogFile::DataSeparator << -c->elongation.c2()
            << LogFile::DataSeparator << c->polarityR.c1() << LogFile::DataSeparator << -c->polarityR.c2() << LogFile::DataSeparator << c->intIntensityR
            << LogFile::DataSeparator << c->polarityG.c1() << LogFile::DataSeparator << -c->polarityG.c2() << LogFile::DataSeparator << c->intIntensityG
            << LogFile::DataSeparator << c->polarityB.c1() << LogFile::DataSeparator << -c->polarityB.c2() << LogFile::DataSeparator << c->intIntensityB
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
