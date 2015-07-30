/* 
 * File:   LogCellNumber.cpp
 * Author: mmpi
 * 
 * Created on October 4, 2013, 6:08 PM
 */

#include "LogCellNumber.h"

#include "TissueState.h"

bool LogCellNumber::startFile(const std::string &Path) {
  if(!LogFile::startFile(Path)) {
    return false;
  }
  // write header line
  _outputFile << "# ";
  _outputFile << "time" << HeaderSeparator;
  _outputFile << "total cell number" << HeaderSeparator;
  _outputFile << std::endl;
  return true;
}


void LogCellNumber::addFrame(const TissueState *s) {
  _outputFile << s->time() << DataSeparator;
  _outputFile << s->numberOfCells();
  _outputFile << std::endl;
}

