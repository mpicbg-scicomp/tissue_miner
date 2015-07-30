/* 
 * File:   LogFileCellNumberChange.cpp
 * Author: mmpi
 * 
 * Created on October 4, 2013, 5:47 PM
 */

#include "LogCellNumberChange.h"
#include "TissueState.h"

bool LogCellNumberChange::startFile(const std::string &Path) {
  if(!LogFile::startFile(Path)) {
    return false;
  }
  // write header line
  _outputFile << "# ";
  _outputFile << "time (mid interval)" << HeaderSeparator;
  _outputFile << "delta t" << HeaderSeparator;
  logQuantityNames();
  _outputFile << std::endl;
  return true;
}


void LogCellNumberChange::addFrameTransition(const TissueState *s, const TissueState *next) {
  _apoptosis = _divisions = _segErrAppearance = _movedIntoMask = _movesOutOfMask = 0;
  for(TissueState::CellConstIterator it = s->beginCellIterator(); it!=s->endCellIterator(); ++it) {
     Cell *c = s->cell(it);
    if(c->duringTransitionAfter==Cell::Divides) {
      ++_divisions;
    } else if(c->duringTransitionAfter==Cell::Apoptosis) {
      ++_apoptosis;
    } else if(c->duringTransitionAfter==Cell::MovesOutOfMask) {
      ++_movesOutOfMask;
    }
  }
  for(TissueState::CellConstIterator it = next->beginCellIterator(); it!=next->endCellIterator(); ++it) {
    Cell *c = next->cell(it);
    if(c->duringTransitionBefore==Cell::SegmentationErrorAppearance) {
      ++_segErrAppearance;
    } else if(c->duringTransitionBefore==Cell::MovedIntoMask) {
      ++_movedIntoMask;
    }
  }
  _outputFile << 0.5*(s->time()+next->time()) << DataSeparator;
  _outputFile << next->time()-s->time() << DataSeparator;
  logQuantities();
  _outputFile << std::endl;
}


void LogCellNumberChange::logQuantities() {
  _outputFile << DataSeparator << _divisions;
  _outputFile << DataSeparator << _apoptosis;
  _outputFile << DataSeparator << _segErrAppearance;
  _outputFile << DataSeparator << _movesOutOfMask;
  _outputFile << DataSeparator << _movedIntoMask;
}
void LogCellNumberChange::logQuantityNames() {
  _outputFile << "divisions" << HeaderSeparator;
  _outputFile << "cell extrusions" << HeaderSeparator;
  _outputFile << "cell appearances due to segmentation errors" << HeaderSeparator;
  _outputFile << "cells moving out of the margin" << HeaderSeparator;
  _outputFile << "cells moving into the margin" << HeaderSeparator;
}
