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

#include "fileIO.h"
#include "LogFile.h"
#include "DeformingTissue.h"

const std::string DeformingTissue::_dbSubFolder("dbTablesFromParser/");


bool DeformingTissue::loadFromSegmentedData(const int MaxFrames) {
  // clean state before
  cleanUp();

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
//    s->markMarginVertices();
    
    // check divisions
    if(!parseDivisionsFile(completeTrackedCellsFileName, completeDaughterCellsFileName, i)) {
      std::cout << "Could not load " << completeDaughterCellsFileName << "!" << std::endl;
      cleanUp();
      return false;
    }
  }

finish:
  if(_frames.size()<2) {
    return false;
  }
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


void DeformingTissue::cleanUp() {
  if(_frames.size()>0) {
    std::cout << "Deleting states..." << std::endl;
    TextProgressBar bar;
    bar.update(0.0);
    for(unsigned int i=0; i<_frames.size() ; ++i) {
      delete _frames[i];
      bar.update((i+1.0)/_frames.size());
    }
    _frames.clear();
    bar.done(true);
  }
}
  
DeformingTissue::~DeformingTissue() {
  cleanUp();
}

std::string DeformingTissue::dbFolderName() const {
  stringstream ss;
  ss << _movie.moviePath << _dbSubFolder;
  return ss.str();
}

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
          << LogFile::HeaderSeparator << "xx component of cell polarity, red channel" << LogFile::HeaderSeparator << "xy component of cell polarity, red channel" << LogFile::HeaderSeparator << "angle integral of intensity, red channel"
          << LogFile::HeaderSeparator << "xx component of cell polarity, green channel" << LogFile::HeaderSeparator << "xy component of cell polarity, green channel" << LogFile::HeaderSeparator << "angle integral of intensity, green channel"
          << LogFile::HeaderSeparator << "xx component of cell polarity, blue channel" << LogFile::HeaderSeparator << "xy component of cell polarity, blue channel" << LogFile::HeaderSeparator << "angle integral of intensity, blue channel"
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
