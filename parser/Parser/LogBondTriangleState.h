/* 
 * File:   LogBondTriangleState.h
 * Author: mmpi
 *
 * Created on October 8, 2013, 7:22 PM
 */

#ifndef LOGBONDTRIANGLESTATE_H
#define	LOGBONDTRIANGLESTATE_H

#include "TissueState.h"
#include "Triangle.h"
#include "LogFile.h"

class LogBondTriangleState : public LogFile {
public:
  LogBondTriangleState(bool logTime);
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
  Nematic2D _wSumQ, _sumQ;
};

#endif	/* LOGBONDTRIANGLESTATE_H */

