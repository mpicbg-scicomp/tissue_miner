/* 
 * File:   DeformingTissue.cpp
 * Author: mmpi
 * 
 * Created on September 25, 2013, 4:48 PM
 */

#include <stdio.h>
#include <stack>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

#include "TextProgressBar.h"
#include "QtRasterImage.h"
#include "QtVectorImage.h"

#include "fileIO.h"
#include "LogCellNumber.h"
#include "LogCellNumberChange.h"
#include "LogDualMarginState.h"
#include "LogDualMarginDeformation.h"
#include "LogTriangleState.h"
#include "LogTriangleDeformation.h"
#include "LogBondTriangleState.h"
#include "DeformingTissue.h"

const std::string DeformingTissue::_dbSubFolder("dbTablesFromParser/");
const std::string DeformingTissue::_dataSubFolder("parserData/");
const std::string DeformingTissue::_imagesSubFolder("parserData/");
const std::string DeformingTissue::_deformationFolder("parserData/");


bool DeformingTissue::loadFromSegmentedData(const int MaxFrames) {
  // clean state before
  cleanUp();

  // new label
  _label = "raw";

  // loop over frames of movie
  int max = ((MaxFrames>=0) && (MaxFrames<_movie.numFrames))?MaxFrames:_movie.numFrames;
  for(int i=0; i<max; ++i) {
    // check, if files exist...
    std::string completeTrackedCellsFileName(_movie.trackedCellsImageFilePath(i));
    if(!fileExists(completeTrackedCellsFileName)) {
      std::cout << "Could not open " << completeTrackedCellsFileName << ". Stopping here." << std::endl;
      goto finish;
      return true;
    }
    std::string completeOriginalFileName(_movie.originalImageFilePath(i));
    if(!fileExists(completeOriginalFileName)) {
      std::cout << "Could not open " << completeOriginalFileName << ". Stopping here." << std::endl;
      goto finish;
      return true;
    }
    std::string completeDaughterCellsFileName(_movie.daughterCellsImageFilePath(i));
    if(!fileExists(completeDaughterCellsFileName)) {
      std::cout << "Could not open " << completeDaughterCellsFileName << ". Stopping here." << std::endl;
      goto finish;
      return true;
    }
    
    // create tissue state
    TissueState *s = new TissueState(i, _movie.frameTime[i]);
    // add
    _frames.push_back(s);
    // load tissue state from image
    std::cout << "Loading " << completeTrackedCellsFileName << "..." << std::endl;
    if(!s->parseFromTrackedCellsFile(completeTrackedCellsFileName, completeOriginalFileName)) {
      std::cout << "Could not load " << completeTrackedCellsFileName << "!" << std::endl;
      cleanUp();
      return false;
    }
    if(!s->checkTopologicalConsistency()) {
      std::cout << "Topological inconsistencies!" << std::endl;
      cleanUp();
      return false;
    }
    s->markMarginVertices();
    
    // check divisions
    if(!parseDivisionsFile(completeTrackedCellsFileName, completeDaughterCellsFileName, i)) {
      std::cout << "Could not load " << completeDaughterCellsFileName << "!" << std::endl;
      cleanUp();
      return false;
    }
  }

finish:
  std::cout << "Removing margin cells..." << std::endl;
  // remove margin cells (they have to be included in the first place in order to correctly interprete the division data from packing analyzer)
  TextProgressBar bar;
  for(unsigned int i=0; i<numberOfFrames(); ++i) {
    TissueState *s = frame(i);
    s->removeMarginCells();
    if(!s->checkTopologicalConsistency()) {
      std::cout << "Topological inconsistencies after removing margin cells!" << std::endl;
      cleanUp();
      return false;
    }
    bar.update(1.0*(i+1)/numberOfFrames());
  }
  bar.done();
  fixCellStatus();
  _label = "rawIncludingCellDivisions";
  return true;
}

int findCellDirection(const PixelFrame &trackedCellsFrame, const Cell *c, const int bondIndex) {
  Pixel p(trackedCellsFrame, trackedCellsFrame.toPosition(c->bonds[bondIndex]->rightVertex->r));
  int foundDirection = -1;
  for(int index=0; index<PixelFrame::NumberOfNeighbors; ++index) {
    if(p.isNeighborOnCanvas(index) && (p.neighbor(index).data()==c->id)) {
      foundDirection = index;
      break;
    }
  }
  return foundDirection;
}

void recursivelyAddToCluster(Cell *c, vector<pair<Cell*,vector<Cell*> > > &mothers, set<Cell*> &visited, set<Cell*> &daughters) {
  if(visited.count(c)==0) {
    // cell is visited
    visited.insert(c);
            
    // add pair, if is mother
    if(c->mother) {
      pair<Cell*,vector<Cell*> > newPair;
      newPair.first = c;
      for(unsigned int i=0; i<c->bonds.size(); ++i) {
        DirectedBond *b = c->bonds[i];
//        if(b->conjBond) {
//          Cell *oc = b->conjBond->cell;
//          if((oc->duringTransitionBefore==Cell::Divided) && (oc->sister==0) && (oc->mother==0)) {
//            newPair.second.push_back(oc);
//          }
//        }
        Vertex *v = b->rightVertex;
        for(unsigned int k=0; k<v->bonds.size(); ++k) {
          DirectedBond *vb = v->bonds[k];
          if(vb!=b) {
            Cell *oc = vb->cell;
            if((oc->duringTransitionBefore==Cell::Divided) && (oc->sister==0) && (oc->mother==0)) {
              newPair.second.push_back(oc);
            }
          }
        }
      }
      mothers.push_back(newPair);
    } else {
      daughters.insert(c);
    }

    // always continue recursion
    for(unsigned int i=0; i<c->bonds.size(); ++i) {
      DirectedBond *b = c->bonds[i];
      Vertex *v = b->rightVertex;
      for(unsigned int k=0; k<v->bonds.size(); ++k) {
        DirectedBond *vb = v->bonds[k];
        if(vb!=b) {
          Cell *oc = vb->cell;
          if((oc->duringTransitionBefore==Cell::Divided) && (oc->sister==0)) {
            recursivelyAddToCluster(oc, mothers, visited, daughters);
          }
        }
      }
    }
  }
}

double findMinimalTotalDistance(vector<pair<Cell*,vector<Cell*> > > &mothers, vector<Cell*> &optimalDaughters, int index, set<Cell*> &occupiedDaughters) {
  if(index==mothers.size()) {
    return 0.0;
  }
  
  double minDistanceSum = -1;
  
  vector<Cell*> &neighbors = mothers[index].second;
  for(unsigned int i=0; i<neighbors.size(); ++i) {
    Cell *c = neighbors[i];
    if(occupiedDaughters.count(c)==0) {
      vector<Cell*> optimalDaughtersParticular(optimalDaughters);
      set<Cell*> occupiedDaughtersParticular(occupiedDaughters);
      optimalDaughtersParticular[index] = c;
      occupiedDaughtersParticular.insert(c);
      double distanceSum = findMinimalTotalDistance(mothers, optimalDaughtersParticular, index+1, occupiedDaughtersParticular);
      if(distanceSum<-0.5) { // not found any combination
        continue;
      }
      distanceSum += (mothers[index].first->r-c->r).norm();
      if((minDistanceSum<-0.5) || (distanceSum<minDistanceSum)) {
        minDistanceSum = distanceSum;
        optimalDaughters.swap(optimalDaughtersParticular);
        occupiedDaughters.swap(occupiedDaughtersParticular);
      }
    }
  }
  
  return minDistanceSum;
}

bool DeformingTissue::parseDivisionsFile(const string TrackedCellsFileName, const string DivisionsFileName, const int frame) {
  if(frame>0) {
    cout << "Parsing " << DivisionsFileName << "..." << endl;
    
    TissueState *sBefore = _frames[frame-1];
    TissueState *sAfter = _frames[frame];

    // get data from image files
    QtRasterImage trackedCellsImage(TrackedCellsFileName.c_str());
    if(!trackedCellsImage.valid()) {
      return false;
    }
    PixelFrame trackedCellsFrame(trackedCellsImage);
    
    QtRasterImage divisionsImage(DivisionsFileName.c_str());
    if(!divisionsImage.valid()) {
      return false;
    }
    PixelFrame divisionsFrame(divisionsImage);
    
    for(TissueState::CellIterator it=sAfter->beginCellIterator(); it!=sAfter->endCellIterator(); ++it) {
      Cell *c = sAfter->cell(it);
      int foundDirection = -1;
      int foundBondIndex = -1;
      for(unsigned int bondIndex=0; bondIndex<c->bonds.size(); ++bondIndex) {
        foundDirection = findCellDirection(trackedCellsFrame, c, bondIndex);
        if(foundDirection>=0) {
          foundBondIndex = bondIndex;
          break;
        }
      }
      if(foundDirection<0) {
        cout << "DeformingTissue::parseDivisionsFile: Cell index " << c->id << " not found!" << endl;
        return false;
      }
      Pixel p(divisionsFrame, divisionsFrame.toPosition(c->bonds[foundBondIndex]->rightVertex->r));
      PixelValue pv = p.neighbor(foundDirection).data();
      if(pv==Pixel::DividingCellValue) {
        c->duringTransitionBefore = Cell::Divided;
        if(sBefore->contains(c->id)) {
          c->mother = c->id;
        }
      } else if(pv!=0) {
        if(pv==Pixel::BondValue) {
          cout << "DeformingTissue::parseDivisionsFile: hit cell boundary at position " << p.positionString() << endl;
          return false;
        } else {
          cout << "DeformingTissue::parseDivisionsFile: hit unknown pixel value " << std::hex << pv << " at position " << p.positionString() << endl;
          return false;
        }
      }
    }

    // guess divisions
    for(TissueState::CellIterator it=sAfter->beginCellIterator(); it!=sAfter->endCellIterator(); ++it) {
      Cell *c = sAfter->cell(it);
      // go through all untouched cells
      if((c->duringTransitionBefore==Cell::Divided) && (c->sister==0)) {
        // create structure for cluster of dividing cells
        vector<pair<Cell*,vector<Cell*> > > mothers;
        set<Cell*> visited, daughters;
        recursivelyAddToCluster(c, mothers, visited, daughters);
        
        // find minimal total distance
        vector<Cell*> optimalDaughters(mothers.size(), NULL);
        set<Cell*> occupiedDaughters;
        if((mothers.size()!=daughters.size()) || (findMinimalTotalDistance(mothers, optimalDaughters, 0, occupiedDaughters)<-0.5)) {
          cout << "DeformingTissue::parseDivisionsFile: Did not find any combination of mother-daughter pairs:" << endl;
          cout << "DeformingTissue::parseDivisionsFile: Mothers: " << endl;
          for(unsigned int i=0; i<mothers.size(); ++i) {
            Vector2D r(mothers[i].first->r);
            r.set(int(r.x()+0.5), trackedCellsImage.height()-1-int(r.y()+0.5));
            cout << "DeformingTissue::parseDivisionsFile:   " << r << " (" << mothers[i].first->id << "): ";
            for(unsigned int k=0; k<mothers[i].second.size(); ++k) {
              cout << mothers[i].second[k]->id << "; ";
            }
            cout << endl;
          }
          cout << "DeformingTissue::parseDivisionsFile: Daughters: ";
          for(set<Cell*>::const_iterator it=daughters.begin(); it!=daughters.end(); ++it) {
            Vector2D r((*it)->r);
            r.set(int(r.x()+0.5), trackedCellsImage.height()-1-int(r.y()+0.5));
            cout << r << " (" << (*it)->id << "); ";
          }
          cout << endl;
          return false;
        }
        
        // set all properties
//        cout << "DeformingTissue::parseDivisionsFile: Cluster" << endl;
        for(unsigned int i=0; i<mothers.size(); ++i) {
          Cell *mAfter = mothers[i].first;
          Cell *dAfter = optimalDaughters[i];
          
          if(sBefore->contains(mAfter->id)) {
            if(dAfter) {
              Cell *mBefore = sBefore->cell(mAfter->id);
    //          cout << "DeformingTissue::parseDivisionsFile: mother before: " << mBefore << "; mother after: " << mAfter << "; daughter after: " << dAfter << endl;

              mBefore->duringTransitionAfter = Cell::Divides;
              mBefore->daughter = dAfter->id;
              mAfter->sister = dAfter->id;
              mAfter->mother = mAfter->id;         
              dAfter->sister = mAfter->id;
              dAfter->mother = mAfter->id;
            } else {
              cout << "DeformingTissue::parseDivisionsFile: Daughter cell not found!";
              return false;
            }
          } else {
            cout << "DeformingTissue::parseDivisionsFile: Mother cell not found in frame before!";
            return false;
          }
        }
      }
    }
  }
  
  return true;
}

bool DeformingTissue::loadCellDivisionData(const std::string FileName, const int skipRows) {
  // check status
  if(0!=_label.compare("raw")) {
    cout << "DeformingTissue::loadCellDivisionData: status is not \"raw\", but \"" << _label << "\"!" << endl;
    return false;
  }
  
  // load file
  const std::string AbsoluteFileName(_movie.movieFilePath(FileName));
  cout << "Loading " << AbsoluteFileName << "..." << endl;
  ifstream file(AbsoluteFileName.c_str());
  if(file.fail()) {
    cout << "DeformingTissue::loadCellDivisionData: Unable to open " << AbsoluteFileName << "!" << endl;
    return false;
  }
  
  long lineCount = 0;
  TextProgressBar bar;

  // skip rows
  for(int r=0;r<skipRows;++r) {
    ++lineCount;
    file.ignore(numeric_limits<streamsize>::max(), '\n');
  }

//  std::stack<pair<CellIndex,CellIndex> > sisterToMotherId; // for divisions, where the daughter with the same id as the mother gets lost
  int lastFrame=-1, lastCellId=-1;
  // loop over the rest
  while(!file.eof()) {
    ++lineCount;
      
    // read data of current line
    int frameId, cellId, sisterId;
    file >> frameId >> cellId >> sisterId;
    if(file.eof()) {
      break;
    }
    // check
    if(file.fail()) {
      cout << "DeformingTissue::loadCellDivisionData: Invalid format in line " << lineCount << "!" << endl;
      cleanUp();
      return false;
    }
    // ignore rest
    file.ignore(numeric_limits<streamsize>::max(), '\n');
    
    // check for new division
    if(sisterId && ((cellId!=lastCellId) || (frameId!=lastFrame))) {
      lastFrame = frameId;
      lastCellId = cellId;
      if((frameId>0) && (frameId<numberOfFrames())) {
        // ignore if sister was there before...
        if(_frames[frameId-1]->contains(sisterId)) {
          std::cout << "DeformingTissue::loadCellDivisionData: Sister cell " << std::hex << sisterId << std::dec << "(" << sisterId << ") of cell " << std::hex << cellId << std::dec << "(" << cellId << ") was already there in frame " << frameId-1 << "! Ignoring..." << std::endl;
          continue;
        }
        // look for cells and do checks
        Cell *cFrameBefore = _frames[frameId-1]->cellCheck(cellId);
        if(!cFrameBefore) {
//          std::cout << "DeformingTissue::loadCellDivisionData: Mother cell " << std::hex << cellId << std::dec << "(" << cellId << ") not found in frame " << frameId-1 << "!" << std::endl;
          continue;
//          cleanUp();
//          return false;
        }
        if(cFrameBefore->duringTransitionAfter) {
          std::cout << "DeformingTissue::loadCellDivisionData: Status of mother cell " << std::hex << cellId << std::dec << " in frame " << frameId << " already set!" << std::endl;
          cleanUp();
          return false;
        }
        Cell *sFrameAfter = _frames[frameId]->cellCheck(sisterId);
        if(!sFrameAfter) {
//          std::cout << "DeformingTissue::loadCellDivisionData: Sister cell " << std::hex << sisterId << std::dec << " not found in frame " << frameId << "!" << std::endl;
          continue;
//          cleanUp();
//          return false;
        }
        if(sFrameAfter->duringTransitionBefore) {
          std::cout << "DeformingTissue::loadCellDivisionData: Status of sister cell " << std::hex << sisterId << std::dec << " in frame " << frameId << " already set!" << std::endl;
          cleanUp();
          return false;
        }
        Cell *cFrameAfter = _frames[frameId]->cellCheck(cellId);
        if(!cFrameAfter) {
//          std::cout << "DeformingTissue::loadCellDivisionData: Cell " << std::hex << cellId << std::dec << " not found in frame " << frameId << "!" << std::endl;
//          std::cout << "DeformingTissue::loadCellDivisionData: Cell " << std::hex << cellId << std::dec << " not found in frame " << frameId << " -> changing id of sister cell to that of mother cell in all following frames." << std::endl;
          // NOTE: have to change id of sister cell to that of mother cell in all following frames; but what happens if there is another division?
          // do it afterwards...
//          sisterToMotherId.push(pair<CellIndex,CellIndex>(sisterId, cellId));
          continue;
//          cleanUp();
//          return false;
        }
        if(cFrameAfter->duringTransitionBefore) {
          std::cout << "DeformingTissue::loadCellDivisionData: Status of cell " << std::hex << cellId << std::dec << "(" << cellId << ") in frame " << frameId << " already set!" << std::endl;
          cleanUp();
          return false;
        }
        // set the data
        cFrameBefore->duringTransitionAfter = Cell::Divides;
        cFrameBefore->daughter = sisterId;
        cFrameAfter->duringTransitionBefore = Cell::Divided;
        cFrameAfter->sister = sisterId;
        cFrameAfter->mother = cellId;         
        sFrameAfter->duringTransitionBefore = Cell::Divided;
        sFrameAfter->sister = cellId;
        sFrameAfter->mother = cellId;
      } else {
        std::cout << "DeformingTissue::loadCellDivisionData: Invalid frame number " << frameId << " in " << AbsoluteFileName << "!" << std::endl;
        cleanUp();
        return false;
      }
    }
    bar.update(1.0*frameId/numberOfFrames());
  }
  bar.update(1.0);

//  while(!sisterToMotherId.empty()) {
//    for(int i=0; i<numberOfFrames(); ++i) {
//      Cell *c = frame(i)->cellCheck(sisterToMotherId.top().first);
//      if(c) {
//        if(frame(i)->contains(sisterToMotherId.top().second)) {
//          std::cout << "DeformingTissue::loadCellDivisionData: Could not replace sister id, since mother id is in frame " << i << "!" << std::endl;
//          cleanUp();
//          return false;
//        }
//        c->id = sisterToMotherId.top().second;
//      }
//    }
//    sisterToMotherId.pop();
//  }
  bar.done(true);
  
  fixCellStatus();
  _label = "rawIncludingCellDivisions";
  return true;
}

void recursivelyMarkAsMovesOutOfMask(Cell *c) {
  c->duringTransitionAfter = Cell::MovesOutOfMask;
  for(unsigned int j=0; j<c->bonds.size(); ++j) {
    if(c->bonds[j]->conjBond) {
      Cell *oc = c->bonds[j]->conjBond->cell;
      if(oc->duringTransitionAfter==Cell::UnclassifiedAfter) {
        recursivelyMarkAsMovesOutOfMask(oc);
      }
    }
  }
}

void recursivelyMarkAsMovedIntoMask(Cell *c) {
  c->duringTransitionBefore = Cell::MovedIntoMask;
  for(unsigned int j=0; j<c->bonds.size(); ++j) {
    if(c->bonds[j]->conjBond) {
      Cell *oc = c->bonds[j]->conjBond->cell;
      if(oc->duringTransitionBefore==Cell::UnclassifiedBefore) {
        recursivelyMarkAsMovedIntoMask(oc);
      }
    }
  }
}


void DeformingTissue::fixCellStatus() {
  cout << "Fixing cell status..." << endl;
  TextProgressBar bar;
  TissueState *before = frame(0);
  for(int f=1; f<numberOfFrames(); ++f) {
    TissueState *after=frame(f);
    // first, check cells which stay
    for(TissueState::CellConstIterator it=before->beginCellIterator(); it!=before->endCellIterator(); ++it) {
      Cell *c = before->cell(it);
      Cell *ac = after->cellCheck(c->id);
      if(ac) {
        if(c->duringTransitionAfter==Cell::UnclassifiedAfter) {
          c->duringTransitionAfter = Cell::Stays;
        }
        if(ac->duringTransitionBefore==Cell::UnclassifiedBefore) {
          ac->duringTransitionBefore = Cell::Stayed;
        }
      }
    }
    // check for cells moving out of mask
    for(TissueState::CellConstIterator it=before->beginCellIterator(); it!=before->endCellIterator(); ++it) {
      Cell *c = before->cell(it);
      if(c->duringTransitionAfter==Cell::UnclassifiedAfter) {
        // is margin cell?
        for(unsigned int j=0; j<c->bonds.size(); ++j) {
          if(!c->bonds[j]->conjBond) {
            // then recursively mark non-staying neighbors
            recursivelyMarkAsMovesOutOfMask(c);
            break;
          }
        }
      }
    }
    // check for cells moving into mask
    for(TissueState::CellConstIterator it=after->beginCellIterator(); it!=after->endCellIterator(); ++it) {
      Cell *c = after->cell(it);
      if(c->duringTransitionBefore==Cell::UnclassifiedBefore) {
        // is margin cell?
        for(unsigned int j=0; j<c->bonds.size(); ++j) {
          if(!c->bonds[j]->conjBond) {
            // then recursively mark non-staying neighbors
            recursivelyMarkAsMovedIntoMask(c);
            break;
          }
        }
      }
    }
    // finally, set the rest of frame before
    for(TissueState::CellConstIterator it=before->beginCellIterator(); it!=before->endCellIterator(); ++it) {
      Cell *c = before->cell(it);
      if(c->duringTransitionAfter==Cell::UnclassifiedAfter) {
        c->duringTransitionAfter = Cell::Apoptosis;
      }
    }
    // finally, set the rest of frame after
    for(TissueState::CellConstIterator it=after->beginCellIterator(); it!=after->endCellIterator(); ++it) {
      Cell *c = after->cell(it);
      if(c->duringTransitionBefore==Cell::UnclassifiedBefore) {
        c->duringTransitionBefore = Cell::SegmentationErrorAppearance;
      }
    }
    
    before = after;
    bar.update((f+1.0)/numberOfFrames());
  }
  bar.done(true);
}

void recursivelyAddToToBeRemoved(std::set<CellIndex> &toBeRemoved, Cell *c) {
  toBeRemoved.insert(c->id);
  for(unsigned int j=0; j<c->bonds.size(); ++j) {
    if(c->bonds[j]->conjBond) {
      Cell *oc = c->bonds[j]->conjBond->cell;
      if(((oc->duringTransitionAfter==Cell::Apoptosis) || (oc->duringTransitionBefore==Cell::SegmentationErrorAppearance) || (oc->duringTransitionAfter==Cell::SegmentationErrorDisappearance)) && (toBeRemoved.count(oc->id)==0)) {
        recursivelyAddToToBeRemoved(toBeRemoved, oc);
      }
    }
  }  
}

void DeformingTissue::removeCellsMiomHelper(std::set<CellIndex> &toBeRemoved) {
  TextProgressBar bar;
  int run = 1;
  do {
    cout << "  Adding daughters..." << endl;
    // add daughter cells
    bar.update(0);
    for(int f=0; f<numberOfFrames(); ++f) {
      TissueState *s=frame(f);
      for(std::set<CellIndex>::iterator it=toBeRemoved.begin(); it!=toBeRemoved.end(); ++it) {
        Cell *c = s->cellCheck(*it);
        if(c && (c->duringTransitionAfter==Cell::Divides)) {
          toBeRemoved.insert(c->daughter);
        }
      }
      bar.update((f+1.0)/numberOfFrames());
    }
    bar.done(true);  
    
    cout << "  Removing " << toBeRemoved.size() << " cells..." << endl;
    // add remove cells
    bar.update(0);
    for(int f=0; f<numberOfFrames(); ++f) {
      TissueState *s=frame(f);
      // first, check cells which stay
      for(std::set<CellIndex>::iterator it=toBeRemoved.begin(); it!=toBeRemoved.end(); ++it) {
        if(s->contains(*it)) {
          s->removeCell(*it);
        }
      }
      bar.update((f+1.0)/numberOfFrames());
    }
    bar.done(true);  
    
    ++run;
    toBeRemoved.clear();
    cout << "  Run " << run << " looking for cells..." << endl;
    bar.update(0);
    for(int f=0; f<numberOfFrames(); ++f) {
      TissueState *s=frame(f);
      for(TissueState::CellConstIterator it=s->beginCellIterator(); it!=s->endCellIterator(); ++it) {
        Cell *c = s->cell(it);
        // count neighbors
        int actualNeighbors = 0;
        for(unsigned int j=0; j<c->bonds.size(); ++j) {
          if(c->bonds[j]->conjBond) {
            ++actualNeighbors;
          }
        }
        // diappearing/appearing cells at the margin
        if((actualNeighbors<c->bonds.size()) && ((c->duringTransitionAfter==Cell::Apoptosis) || (c->duringTransitionBefore==Cell::SegmentationErrorAppearance) || (c->duringTransitionAfter==Cell::SegmentationErrorDisappearance)) && (toBeRemoved.count(c->id)==0)) {
          recursivelyAddToToBeRemoved(toBeRemoved, c);
        }
        // cells with at most 1 neighbor
        if(actualNeighbors<2) {
          toBeRemoved.insert(c->id);          
        }
      }
      bar.update((f+1.0)/numberOfFrames());
    }
    bar.done(true);  
  } while(!toBeRemoved.empty());
  cout << "Done." << endl;
  cout << "Updating margin..." << endl;
  bar.update(0);
  for(int f=0; f<numberOfFrames(); ++f) {
    frame(f)->markMarginVertices();
    bar.update((f+1.0)/numberOfFrames());
  }
  bar.done(true);  
  
}

void DeformingTissue::removeCellsMiom() {
  cout << "Removing cells moving into or out of the margin..." << endl;
  std::set<CellIndex> toBeRemoved;

  cout << "  First run looking for cells..." << endl;
  TextProgressBar bar;
  // directly get cells moving into or out of the tissue
  for(int f=0; f<numberOfFrames(); ++f) {
    TissueState *s=frame(f);
    // first, check cells which stay
    for(TissueState::CellConstIterator it=s->beginCellIterator(); it!=s->endCellIterator(); ++it) {
      Cell *c = s->cell(it);
      if((c->duringTransitionAfter==Cell::MovesOutOfMask) || (c->duringTransitionBefore==Cell::MovedIntoMask)) {
        toBeRemoved.insert(c->id);
      }
    }
    bar.update((f+1.0)/numberOfFrames());
  }
  bar.done(true);  
  
  removeCellsMiomHelper(toBeRemoved);
  _label = "removedMiom";
}

bool DeformingTissue::applyMaskToFinalState(const std::string MaskTag) {
  cout << "Applying \"" << MaskTag << "\" mask..." << endl;
  stringstream maskImageFileName;
  maskImageFileName << "masks/" << MaskTag << ".png";
  
  // load mask image
  QtRasterImage maskImage(_movie.segmentedFrameFilePath(_frames.size()-1, maskImageFileName.str()).c_str());
  if(!maskImage.valid()) {
    return false;
  }
  
  // check cells to be removed in last frame
  TissueState *lastFrame = frame(numberOfFrames()-1);
  std::set<CellIndex> cellsToBeRemoved;
  for(TissueState::CellConstIterator it=lastFrame->beginCellIterator(); it!=lastFrame->endCellIterator(); ++it) {
    Cell *c = lastFrame->cell(it);
    Vector2D p(c->r);
    if(!maskImage.pixel(p.x(), p.y(), RedChannel)) {
      cellsToBeRemoved.insert(c->id);
    }
  }
  
  // iteratively remove cells from preceding frames
  removeCellsMiomHelper(cellsToBeRemoved);
  
  // remove cells moving into or out of the margin
  removeCellsMiom();
  
  _label = MaskTag;  
  return true;
}

bool DeformingTissue::writeDeformationData() {
  // create folder if it doesn't exist
  if(!folderExists(textDataFolderName())) {
    if(!recursivelyCreateFolder(textDataFolderName())) {
      std::cout << "Could not create " << textDataFolderName() << "!" << std::endl;
      return false;
    }
  }
  
  TextProgressBar bar;
  
  // write cell number balance data files
  std::cout << "Writing cell number data..." << std::endl;
  LogCellNumber cellNumber;
  if(!cellNumber.startFile(textDataFileName("cellNumber.dat"))) {
    return false;
  }
  bar.update(0);
  for(unsigned int i=0; i<numberOfFrames(); ++i) {
    cellNumber.addFrame(frame(i));
    bar.update((i+1.0)/numberOfFrames());
  }
  if(!cellNumber.endFile()) {
    return false;
  }
  bar.done(true);

  std::cout << "Writing cell number change data..." << std::endl;
  LogCellNumberChange cellNumberBalance;
  if(!cellNumberBalance.startFile(textDataFileName("cellNumberChange.dat"))) {
    return false;
  }
  bar.update(0);
  for(unsigned int i=0; i<numberOfFrames()-1; ++i) {
    cellNumberBalance.addFrameTransition(frame(i), frame(i+1));
    bar.update((i+1.0)/(numberOfFrames()-1));
  }
  if(!cellNumberBalance.endFile()) {
    return false;
  }
  bar.done(true);

  // create intermediate states
  if(_statesBeforeDivision.size() || _statesAfterDeformation.size()) {
    std::cout << "DeformingTissue::writeDeformationData: intermediate states already present!" << std::endl;
    throw std::exception();
  }
  std::cout << "Creating intermediate states..." << std::endl;
  bar.update(0);
  for(unsigned int i=0; i<numberOfFrames()-1; ++i) {
    _statesBeforeDivision.push_back(frame(i+1)->createCopyFusingDivisions(*frame(i)));
    if(!_statesBeforeDivision[i]->checkTopologicalConsistency()) {
      return false;
    }
    _statesBeforeDivision[i]->markMarginVertices();
    _statesAfterDeformation.push_back(frame(i)->createCopyAndMoveCellPositionsTo(*_statesBeforeDivision[i]));
    if(!_statesAfterDeformation[i]->checkTopologicalConsistency()) {
      return false;
    }
    _statesAfterDeformation[i]->markMarginVertices();
    bar.update((i+1.0)/(numberOfFrames()-1));
  }
  bar.done(true);

  // creating triangles
  std::cout << "Creating triangles..." << std::endl;
  bar.update(0);
  for(unsigned int i=0; i<numberOfFrames(); ++i) {
    frame(i)->createTriangles();
    if(i<_statesBeforeDivision.size()) {
      _statesBeforeDivision[i]->createTriangles();
    }
    bar.update((i+1.0)/(numberOfFrames()));
  }
  bar.done(true);

  // write data files
  std::cout << "Writing dual margin state data..." << std::endl;
  LogDualMarginState dualMarginState;
  if(!dualMarginState.startFile(textDataFileName("dualMarginState.dat"))) {
    return false;
  }
  bar.update(0);
  for(unsigned int i=0; i<numberOfFrames(); ++i) {
    dualMarginState.addState(frame(i));
    bar.update((i+1.0)/numberOfFrames());
  }
  if(!dualMarginState.endFile()) {
    return false;
  }
  bar.done(true);

  std::cout << "Writing dual margin deformation data..." << std::endl;
  LogDualMarginDeformation dualMarginDeformation;
  if(!dualMarginDeformation.startFile(textDataFileName("dualMarginDeformation.dat"))) {
    return false;
  }
  bar.update(0);
  for(unsigned int i=0; i<numberOfFrames()-1; ++i) {
    dualMarginDeformation.addDeformation(frame(i), frame(i+1));
    bar.update((i+1.0)/(numberOfFrames()-1));
  }
  if(!dualMarginDeformation.endFile()) {
    return false;
  }
  bar.done(true);

  std::cout << "Writing average triangle state data..." << std::endl;
  LogTriangleState avgState(true);
  if(!avgState.startFile(textDataFileName("triangleState.dat"))) {
    return false;
  }
  bar.update(0);
  for(unsigned int i=0; i<numberOfFrames(); ++i) {
    avgState.addState(frame(i));
    bar.update((i+1.0)/numberOfFrames());
  }
  if(!avgState.endFile()) {
    return false;
  }
  bar.done(true);

  std::cout << "Writing average bond triangle state data..." << std::endl;
  LogBondTriangleState avgBondTriangleState(true);
  if(!avgBondTriangleState.startFile(textDataFileName("bondTriangleState.dat"))) {
    return false;
  }
  bar.update(0);
  for(unsigned int i=0; i<numberOfFrames(); ++i) {
    avgBondTriangleState.addState(frame(i));
    bar.update((i+1.0)/numberOfFrames());
  }
  if(!avgBondTriangleState.endFile()) {
    return false;
  }
  bar.done(true);

  std::cout << "Writing average triangle state before division data..." << std::endl;
  LogTriangleState avgStateBeforeDivision(false);
  if(!avgStateBeforeDivision.startFile(textDataFileName("triangleStateBeforeDivision.dat"))) {
    return false;
  }
  bar.update(0);
  for(unsigned int i=0; i<_statesBeforeDivision.size(); ++i) {
    avgStateBeforeDivision.addState(_statesBeforeDivision[i]);
    bar.update((i+1.0)/_statesBeforeDivision.size());
  }
  if(!avgStateBeforeDivision.endFile()) {
    return false;
  }
  bar.done(true);
  
  std::cout << "Writing average triangle state after deformation data..." << std::endl;
  LogTriangleState avgStateAfterDeformation(false);
  if(!avgStateAfterDeformation.startFile(textDataFileName("triangleStateAfterDeformation.dat"))) {
    return false;
  }
  bar.update(0);
  for(unsigned int i=0; i<_statesAfterDeformation.size(); ++i) {
    avgStateAfterDeformation.addState(_statesAfterDeformation[i]);
    bar.update((i+1.0)/_statesAfterDeformation.size());
  }
  if(!avgStateAfterDeformation.endFile()) {
    return false;
  }
  bar.done(true);
  
  std::cout << "Writing average triangle deformation data..." << std::endl;
  LogTriangleDeformation avgDeformation;
  if(!avgDeformation.startFile(textDataFileName("triangleDeformation.dat"))) {
    return false;
  }
  bar.update(0);
  for(unsigned int i=0; i<numberOfFrames()-1; ++i) {
    avgDeformation.addDeformation(frame(i), frame(i+1));
    bar.update((i+1.0)/(numberOfFrames()-1));
  }
  if(!avgDeformation.endFile()) {
    return false;
  }
  bar.done(true);
  
  return true;
}


bool DeformingTissue::drawFrames(const std::string Tag, void (*drawFrame)(QtImage &img, const IsotropicTransformation &Ref2Pixel, const TissueState& s), const bool smaller) const {
  if(!folderExists(imageFolderName(Tag))) {
    recursivelyCreateFolder(imageFolderName(Tag));
  }
  const int SmallerWidth = 1024;
  const int SmallerHeight = 768;
  for(int i=0; i<numberOfFrames(); ++i) {
//    QtRasterImage img(_movie.width, _movie.height);
//    img.reset(Black);
    QtRasterImage img(_movie.originalImageFilePath(frame(i)->frameNumber()).c_str());
    if(!img.valid()) {
      return false;
    }
    drawFrame(img, Identity, *frame(i));
    if(smaller) {
      img.scaled(SmallerWidth, SmallerHeight).save(imageFileName(Tag, ".png", frame(i)->frameNumber()).c_str());
    } else {
      img.save(imageFileName(Tag, ".png", frame(i)->frameNumber()).c_str());
    }
  }
  return true;
}

bool DeformingTissue::drawFramesTcVec(const std::string Tag, void (*drawFrame)(QtImage &img, const IsotropicTransformation &Ref2Pixel, const TissueState& s)) const {
  if(!folderExists(imageFolderName(Tag))) {
    recursivelyCreateFolder(imageFolderName(Tag));
  }
  for(int i=0; i<numberOfFrames(); ++i) {
    QtRasterImage img(_movie.trackedCellsImageFilePath(frame(i)->frameNumber()).c_str());
    if(!img.valid()) {
      return false;
    }
    QtVectorImage vecImg(img);
    if(!vecImg.start(imageFileName(Tag, ".pdf", frame(i)->frameNumber()).c_str())) {
      return false;
    }
    vecImg.drawImage(img);
    drawFrame(vecImg, Identity, *frame(i));
    vecImg.done();
  }
  return true;
}

bool DeformingTissue::drawStatesBeforeDivisionDebug() const {
  std::string Tag("statesBeforeDivisionDebug");
  if(!folderExists(imageFolderName(Tag))) {
    recursivelyCreateFolder(imageFolderName(Tag));
  }
  IsotropicTransformation trafo(IsotropicTransformation(Vector2D(0.5, -0.5)));
  for(int i=0; i<_statesBeforeDivision.size(); ++i) {
    QtRasterImage img(_movie.trackedCellsImageFilePath(i+1).c_str());
    if(!img.valid()) {
      return false;
    }
    QtVectorImage vecImg(img);
    if(!vecImg.start(imageFileName(Tag, ".pdf", i).c_str())) {
      return false;
    }
    vecImg.drawImage(img);
    _statesBeforeDivision[i]->drawBonds(vecImg, trafo, Black, 0.5);
    _statesBeforeDivision[i]->drawVertices(vecImg, trafo, Red, 0.5);
    _statesBeforeDivision[i]->drawCells(vecImg, trafo, White, 3);
    _statesBeforeDivision[i]->drawCellsTransitionAfter(vecImg, trafo, 2);
    vecImg.done();
  }
  return true;
}


void DeformingTissue::cleanUp() {
  if(_frames.size()>0) {
    std::cout << "Deleting states..." << std::endl;
    TextProgressBar bar;
    bar.update(0.0);
    for(unsigned int i=0; i<_frames.size() ; ++i) {
      delete _frames[i];
      if(i<_statesBeforeDivision.size()) {
        delete _statesBeforeDivision[i];      
      }
      if(i<_statesAfterDeformation.size()) {
        delete _statesAfterDeformation[i];      
      }
      bar.update((i+1.0)/_frames.size());
    }
    _frames.clear();
    bar.done(true);
  }
  _label = "void";
}
  
DeformingTissue::~DeformingTissue() {
  cleanUp();
}

std::string DeformingTissue::textDataFolderName() const {
  stringstream ss;
  ss << _movie.moviePath << _deformationFolder << _label << "/";
  return ss.str();
}

std::string DeformingTissue::textDataFileName(const std::string &FileName) const {
  stringstream ss;
  ss << textDataFolderName() << FileName;
  return ss.str();
}

std::string DeformingTissue::imageFolderName(const std::string &Tag) const {
  stringstream ss;
  ss << _movie.moviePath << _imagesSubFolder << _label << "/" << Tag << "/";
  return ss.str();
}

std::string DeformingTissue::imageFileName(const std::string &Tag, const std::string &Extension, const int frame) const {
  stringstream ss;
  ss << imageFolderName(Tag) << "frame_" << std::setw(4) << std::setfill('0') << frame << Extension;
  return ss.str();
}
  
std::string DeformingTissue::dataFolderName() const {
  stringstream ss;
  ss << _movie.moviePath << _dataSubFolder << _label << "/data/";
  return ss.str();
}

std::string DeformingTissue::dataFileName(const int frame) const {
  stringstream ss;
  ss << dataFolderName() << "data" << std::setw(4) << std::setfill('0') << frame << ".nc";
  return ss.str();
}

std::string DeformingTissue::dbFolderName() const {
  stringstream ss;
  ss << _movie.moviePath << _dbSubFolder;
  return ss.str();
}

#ifdef USE_NETCDF

DeformingTissue::DeformingTissue(const MovieData &movie, const std::string label, const int maxFrames) : _movie(movie), _label(label) {
  load(label, maxFrames);
}

bool DeformingTissue::load(const std::string label, const int maxFrames) {
  // set label
  _label = label;
  
  // check out how many files are there...
  int frames = -1;
  do {
    ++frames;
  } while(fileExists(dataFileName(frames)));
  
  // max number of frames
  if((maxFrames>=0) && (frames>maxFrames)) {
    frames = maxFrames;
  }
  
  // load files
  std::cout << "Loading " << frames << " data files from " << dataFolderName() << "..." << std::endl;
  TextProgressBar bar;
  bar.update(0);
  for(int frame=0; frame<frames; ++frame) {
//     cout << "frame " << frame << endl;
    // create tissue state
    TissueState *s = new TissueState;
    // add
    _frames.push_back(s);
    // load tissue state from data file
    if(!s->load(dataFileName(frame))) {
      std::cout << "Could not load " << dataFileName(frame) << "!" << std::endl;
      cleanUp();
      return false;
    }
    bar.update((frame+1.0)/frames);
  }
  bar.done(true);
  return true;
}

bool DeformingTissue::save() const {
  std::cout << "Saving data files to " << dataFolderName() << "..." << std::endl;
  
  // create folder if necessary
  if(folderExists(dataFolderName())) {
    if(!removeFolderContent(dataFolderName())) {
      std::cout << "Could not remove old content of " << dataFolderName() << "!" << std::endl;
      return false;
    }
  }
  if(!folderExists(dataFolderName())) {
    if(!recursivelyCreateFolder(dataFolderName())) {
      std::cout << "Could not create " << dataFolderName() << "!" << std::endl;
      return false;
    }
  }
  
  TextProgressBar bar;
  bar.update(0);
  for(unsigned int i=0; i<_frames.size(); ++i) {
    if(!_frames[i]->save(dataFileName(i))) {
      return false;
    }
    bar.update((i+1.0)/_frames.size());
  }
  bar.done(true);
  return true;
}

#endif  


bool DeformingTissue::exportToDbTables() const {
  std::cout << "Exporting db tables to " << dbFolderName() << "..." << std::endl;
  
  // create folder if necessary
  if(folderExists(dbFolderName())) {
    if(!removeFolderContent(dbFolderName())) {
      std::cout << "Could not remove old content of " << dbFolderName() << "!" << std::endl;
      return false;
    }
  }
  if(!folderExists(dbFolderName())) {
    if(!recursivelyCreateFolder(dbFolderName())) {
      std::cout << "Could not create " << dbFolderName() << "!" << std::endl;
      return false;
    }
  }
  
  // create files
  stringstream fileName("");
  fileName << dbFolderName() << "frame.dat";
  ofstream framesFile(fileName.str().c_str());
  framesFile << "# " << "frame number" << LogFile::HeaderSeparator << "frame time [hAPF]" << LogFile::HeaderSeparator << "\n";
  
  fileName.str("");
  fileName << dbFolderName() << "vertex_in_frame.dat";
  ofstream verticesFile(fileName.str().c_str());
  verticesFile << "# " << "frame number" << LogFile::HeaderSeparator << "vertex_in_frame index" << LogFile::HeaderSeparator << "vertex position x (pixels)" << LogFile::HeaderSeparator << " vertex position y (pixels, zero is topmost pixel row of image)" << LogFile::HeaderSeparator << "\n";
  
  fileName.str("");
  fileName << dbFolderName() << "cell_in_frame.dat";
  ofstream cellsFile(fileName.str().c_str());
  cellsFile << "# " << "frame number" << LogFile::HeaderSeparator << "tracked cell id"
          << LogFile::HeaderSeparator << "during transition before" << LogFile::HeaderSeparator << "during transition after" << LogFile::HeaderSeparator << "daughter cell id"
          << LogFile::HeaderSeparator << "center x (pixel)" << LogFile::HeaderSeparator << "center y (pixel, zero is topmost pixel row of image)" << LogFile::HeaderSeparator << "area"
          << LogFile::HeaderSeparator << "xx component of cell elongation" << LogFile::HeaderSeparator << "xy component of cell elongation (y axis point upwards, here!)"
          << LogFile::HeaderSeparator << "xx component of polarity, red channel" << LogFile::HeaderSeparator << "xy component of cell elongation (y axis point upwards, here!), red channel" << LogFile::HeaderSeparator << "angle integral of intensity, red channel"
          << LogFile::HeaderSeparator << "xx component of polarity, green channel" << LogFile::HeaderSeparator << "xy component of cell elongation (y axis point upwards, here!), green channel" << LogFile::HeaderSeparator << "angle integral of intensity, green channel"
          << LogFile::HeaderSeparator << "xx component of polarity, blue channel" << LogFile::HeaderSeparator << "xy component of cell elongation (y axis point upwards, here!), blue channel" << LogFile::HeaderSeparator << "angle integral of intensity, blue channel"
          << "\n";

  fileName.str("");
  fileName << dbFolderName() << "ignoredCellsAtMargin_in_frame.dat";
  ofstream ignoredCellsFile(fileName.str().c_str());
  ignoredCellsFile << "# " << "frame number" << LogFile::HeaderSeparator << "tracked cell id" << "\n";
  
  fileName.str("");
  fileName << dbFolderName() << "undirectedBond_in_frame.dat";
  ofstream undirectedBondsFile(fileName.str().c_str());
  undirectedBondsFile << "# " << "frame number" << LogFile::HeaderSeparator << "undirectedBond_in_frame index" << LogFile::HeaderSeparator << "bond length (curved, in pixels)" << LogFile::HeaderSeparator << "\n";
  
  fileName.str("");
  fileName << dbFolderName() << "directedBond_in_frame.dat";
  ofstream directedBondsFile(fileName.str().c_str());
  directedBondsFile << "# " << "frame number" << LogFile::HeaderSeparator << "directedBond_in_frame index" << LogFile::HeaderSeparator << "conjugated directedBond_in_frame index" << LogFile::HeaderSeparator << "undirectedBond_in_frame index" << LogFile::HeaderSeparator << "tracked cell id" << LogFile::HeaderSeparator << "vertex_in_frame index" << LogFile::HeaderSeparator << "index of directedBond_in_frame left of bond as seen from cell" << LogFile::HeaderSeparator << "\n";

  unsigned long lastVid=0, lastDbid=0, lastUbid=0;
 
  TextProgressBar bar;
  bar.update(0);
  for(unsigned int i=0; i<_frames.size(); ++i) {
    QtRasterImage img(_movie.trackedCellsImageFilePath(i).c_str());
    if(!_frames[i]->exportToDbTables(img.height(), framesFile, verticesFile, cellsFile, ignoredCellsFile, undirectedBondsFile, directedBondsFile, lastVid, lastDbid, lastUbid)) {
      return false;
    }
    bar.update((i+1.0)/_frames.size());
  }
  
  framesFile.close();
  verticesFile.close();
  cellsFile.close();
  ignoredCellsFile.close();
  undirectedBondsFile.close();
  directedBondsFile.close();

  fileName.str("");
  fileName << dbFolderName() << "parserVersion.dat";
  ofstream versionFile(fileName.str().c_str());
  versionFile << PROGRAMNAME << ", Version: " << VERSION << "\n";
  versionFile.close();

  bar.done(true);
  return true;
}

bool DeformingTissue::exportAllCellIds() const {
  // find all cell ids
  TextProgressBar bar;
  std::cout << "Find all cell ids..." << std::endl;
  bar.update(0);
  set<CellIndex> allCellIds;
  for(unsigned int i=0; i<_frames.size(); ++i) {
    TissueState *s = _frames[i];
    for(TissueState::CellConstIterator it=s->beginCellIterator(); it!=s->endCellIterator(); ++it) {
      allCellIds.insert(s->cell(it)->id);
    }
    bar.update((i+1.0)/_frames.size());
  }  
  bar.done(true);

  
  // create filename
  stringstream fileName("");
  fileName << textDataFolderName() << "allCellIds.dat";
  std::cout << "Exporting all cell ids to " << fileName.str() << "..." << std::endl;

  // create folder if necessary
  if(!folderExists(textDataFolderName())) {
    if(!recursivelyCreateFolder(textDataFolderName())) {
      std::cout << "Could not create " << textDataFolderName() << "!" << std::endl;
      return false;
    }
  }

  // create file
  ofstream cellFile(fileName.str().c_str());
  cellFile << "# " << "cell id" << "\n";
  cellFile << Cell::VoidCellId << std::endl;
  bar.update(0);
  double sum=0.0, offset=1.0/allCellIds.size();
  for(set<CellIndex>::const_iterator it=allCellIds.begin(); it!=allCellIds.end(); ++it) {
    cellFile << *it << std::endl;
    sum += offset;
    bar.update(sum);
  }  
  bar.done(true);
  cellFile.close();
  return true;
}
