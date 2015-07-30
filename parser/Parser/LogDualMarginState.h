/* 
 * File:   LogMarginState.h
 * Author: mmpi
 *
 * Created on October 1, 2013, 8:53 PM
 */

#ifndef LOGDUALMARGINSTATE_H
#define	LOGDUALMARGINSTATE_H

#include "TissueState.h"
#include "Triangle.h"
#include "LogFile.h"

class LogDualMarginState : public LogFile {
public:
  virtual bool startFile(const std::string &Path);
  void addState(const TissueState *s);
};

#endif	/* LOGMARGINSTATE_H */

