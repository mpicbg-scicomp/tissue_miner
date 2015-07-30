/* 
 * File:   TriangleStateAverage.h
 * Author: mmpi
 *
 * Created on 30. September 2013, 18:52
 */

#ifndef LOGTRIANGLESTATE_H
#define	LOGTRIANGLESTATE_H

#include "TissueState.h"
#include "Triangle.h"
#include "LogFile.h"

class LogTriangleState : public LogFile {
public:
  LogTriangleState(bool logTime);
  virtual bool startFile(const std::string &Path);
  void addState(const TissueState *s);
  
private:
  void logQuantityNames();
  void prepareQuantities();
  void addTriangle(const Triangle *t);
  void logQuantities();

  bool _logTime;
  int _number;
  double _sumArea;
  Nematic2D _wSumQ;
};


#endif	/* TRIANGLESTATEAVERAGE_H */

