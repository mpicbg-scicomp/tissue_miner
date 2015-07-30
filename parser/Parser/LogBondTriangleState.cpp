/* 
 * File:   LogBondTriangleState.cpp
 * Author: mmpi
 * 
 * Created on October 8, 2013, 7:22 PM
 */

#include "LogFile.h"
#include "LogBondTriangleState.h"

LogBondTriangleState::LogBondTriangleState(bool logTime) : _logTime(logTime) {}

bool LogBondTriangleState::startFile(const std::string &Path) {
  if(!LogFile::startFile(Path)) {
    return false;
  }
  // write header line
  _outputFile << "# ";
  if(_logTime) {
    _outputFile << "time" << HeaderSeparator;
  }
  logQuantityNames();
  _outputFile << std::endl;
  return true;
}
void LogBondTriangleState::addState(const TissueState *s) {
  prepareQuantities();
  for(TissueState::BondConstIterator it=s->beginBondIterator(); it!=s->endBondIterator(); ++it) {
    addTriangle(s->bond(it)->bondTriangle);
  }
  if(_logTime) {
    _outputFile << s->time() << DataSeparator;
  }
  logQuantities();
  _outputFile << std::endl;
}

void LogBondTriangleState::addTriangle(const Triangle *t) {
  ++_number;
  _sumArea += t->area;
  _wSumQ += t->area * t->Q;
  _sumQ += t->Q;
}
void LogBondTriangleState::prepareQuantities() {
  _number = 0;
  _sumArea = 0;
  _wSumQ.set(0,0);
  _sumQ = _wSumQ;
}
void LogBondTriangleState::logQuantities() {
  _outputFile << _number;
  _outputFile << DataSeparator << _sumArea;
  _wSumQ /= _sumArea;
  _outputFile << DataSeparator << _wSumQ.c1() << DataSeparator << _wSumQ.c2();
  _sumQ /= _number;
  _outputFile << DataSeparator << _sumQ.c1() << DataSeparator << _sumQ.c2();
}
void LogBondTriangleState::logQuantityNames() {
  _outputFile << "total number" << HeaderSeparator;
  _outputFile << "sum area" << HeaderSeparator;
  _outputFile << "<Q_1>_t" << HeaderSeparator << "<Q_2>_t" << HeaderSeparator;
  _outputFile << "equally weighted average Q_1" << HeaderSeparator << "equally weighted average Q_2" << HeaderSeparator;
}
