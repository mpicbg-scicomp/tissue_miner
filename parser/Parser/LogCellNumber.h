/* 
 * File:   LogCellNumber.h
 * Author: mmpi
 *
 * Created on October 4, 2013, 6:08 PM
 */

#ifndef LOGCELLNUMBER_H
#define	LOGCELLNUMBER_H

#include "TissueState.h"
#include "LogFile.h"

class LogCellNumber : public LogFile {
public:
  virtual bool startFile(const std::string &Path);
  void addFrame(const TissueState *s);
};

#endif	/* LOGCELLNUMBER_H */

