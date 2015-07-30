/* 
 * File:   LogFileCellNumberChange.h
 * Author: mmpi
 *
 * Created on October 4, 2013, 5:47 PM
 */

#ifndef LOGCELLNUMBERCHANGE_H
#define	LOGCELLNUMBERCHANGE_H

#include "TissueState.h"
#include "LogFile.h"

class LogCellNumberChange : public LogFile {
public:
  virtual bool startFile(const std::string &Path);
  void addFrameTransition(const TissueState *s, const TissueState *next);

private:
  void logQuantityNames();
  void logQuantities();

  int _divisions, _apoptosis, _movesOutOfMask;
  int _segErrAppearance, _movedIntoMask;
};

#endif	/* LOGFILECELLNUMBERCHANGE_H */

