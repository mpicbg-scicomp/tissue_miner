/* 
 * File:   LogFile.h
 * Author: mmpi
 *
 * Created on October 1, 2013, 7:07 PM
 */

#ifndef LOGFILE_H
#define	LOGFILE_H

#include <fstream>
#include <string>

class LogFile {
public:
  virtual bool startFile(const std::string &Path);
  virtual bool endFile();
  
  const static std::string HeaderSeparator;
  const static std::string DataSeparator;
  
protected:
  std::fstream _outputFile;
};

#endif	/* LOGFILE_H */

