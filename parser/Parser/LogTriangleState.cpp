#include "LogFile.h"
#include "LogTriangleState.h"

LogTriangleState::LogTriangleState(bool logTime) : _logTime(logTime) {}

bool LogTriangleState::startFile(const std::string &Path) {
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
void LogTriangleState::addState(const TissueState *s) {
  prepareQuantities();
  for(unsigned int j=0; j<s->numberOfVertices(); ++j) {
    Vertex *v = s->vertex(j);
    for(unsigned int k=0; k<v->triangles.size(); ++k) {
      addTriangle(v->triangles[k]);
    }
  }
  if(_logTime) {
    _outputFile << s->time() << DataSeparator;
  }
  logQuantities();
  _outputFile << std::endl;
}

void LogTriangleState::addTriangle(const Triangle *t) {
  ++_number;
  _sumArea += t->area;
  _wSumQ += t->area * t->Q;
}
void LogTriangleState::prepareQuantities() {
  _number = 0;
  _sumArea = 0;
  _wSumQ.set(0,0);
}
void LogTriangleState::logQuantities() {
  _outputFile << _number;
  _outputFile << DataSeparator << _sumArea;
  _wSumQ /= _sumArea;
  _outputFile << DataSeparator << _wSumQ.c1() << DataSeparator << _wSumQ.c2();
}
void LogTriangleState::logQuantityNames() {
  _outputFile << "total number" << HeaderSeparator;
  _outputFile << "sum area" << HeaderSeparator;
  _outputFile << "<Q_1>_t" << HeaderSeparator << "<Q_2>_t" << HeaderSeparator;
}
