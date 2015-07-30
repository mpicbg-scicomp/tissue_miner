/* 
 * File:   LogFile.cpp
 * Author: mmpi
 * 
 * Created on October 1, 2013, 7:07 PM
 */

#include <iostream>
#include <fstream>
#include "LogFile.h"

const std::string LogFile::HeaderSeparator(";\t");
const std::string LogFile::DataSeparator("\t");

bool LogFile::startFile(const std::string &Path) {
  if(_outputFile.is_open()) {
    std::cout << "LogFile::startFile: file was already open!" << std::endl;
    throw std::exception();
  }
  _outputFile.open(Path.c_str(), std::ios_base::out);
  if(!_outputFile.is_open()) {
    std::cout << "LogFile::startFile: Could not open " << Path << " for writing!" << std::endl;
    return false;
  }
  return true;
}

bool LogFile::endFile() {
  _outputFile.close();
  if(_outputFile.fail()) {
    std::cout << "LogFile::endFile: Could not close file!" << std::endl;
    return false;
  }
  return true;
}
