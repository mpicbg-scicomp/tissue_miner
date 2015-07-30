/* 
 * File:   LogDualMarginDeformation.h
 * Author: mmpi
 *
 * Created on October 1, 2013, 9:09 PM
 */

#ifndef LOGDUALMARGINDEFORMATION_H
#define	LOGDUALMARGINDEFORMATION_H

#include "TissueState.h"
#include "Triangle.h"
#include "LogFile.h"

class LogDualMarginDeformation : public LogFile {
public:
  virtual bool startFile(const std::string &Path);
  void addDeformation(const TissueState *s, const TissueState *next);
};

#endif	/* LOGDUALMARGINDEFORMATION_H */

