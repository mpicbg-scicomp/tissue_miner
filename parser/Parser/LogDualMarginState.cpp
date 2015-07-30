/* 
 * File:   LogMarginState.cpp
 * Author: mmpi
 * 
 * Created on October 1, 2013, 8:53 PM
 */

#include "LogDualMarginState.h"

bool LogDualMarginState::startFile(const std::string &Path) {
  if(!LogFile::startFile(Path)) {
    return false;
  }
  // write header line
  _outputFile << "# ";
  _outputFile << "time" << HeaderSeparator;
  _outputFile << "area" << HeaderSeparator;
  _outputFile << std::endl;
  return true;
}

void LogDualMarginState::addState(const TissueState *s) {
  _outputFile << s->time();
  _outputFile << DataSeparator << s->totalAreaOfDualMargin();
  _outputFile << std::endl;
}
