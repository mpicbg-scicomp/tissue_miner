/* 
 * File:   LogDualMarginDeformation.cpp
 * Author: mmpi
 * 
 * Created on October 1, 2013, 9:09 PM
 */

#include "LogDualMarginDeformation.h"

bool LogDualMarginDeformation::startFile(const std::string &Path) {
  if(!LogFile::startFile(Path)) {
    return false;
  }
  // write header line
  _outputFile << "# ";
  _outputFile << "time (mid interval)" << HeaderSeparator;
  _outputFile << "delta t" << HeaderSeparator;
  _outputFile << "Delta u_{xx}" << HeaderSeparator;
  _outputFile << "Delta u_{xy}" << HeaderSeparator;
  _outputFile << "Delta u_{yx}" << HeaderSeparator;
  _outputFile << "Delta u_{yy}" << HeaderSeparator;
  _outputFile << std::endl;
  return true;
}

void LogDualMarginDeformation::addDeformation(const TissueState *s, const TissueState *next) {
  _outputFile << 0.5*(s->time()+next->time());
  _outputFile << DataSeparator << next->time()-s->time();
  Matrix2x2 du(s->totalDeformationOfDualMargin());
  _outputFile << DataSeparator << du.m11();
  _outputFile << DataSeparator << du.m12();
  _outputFile << DataSeparator << du.m21();
  _outputFile << DataSeparator << du.m22();
  _outputFile << std::endl;
}