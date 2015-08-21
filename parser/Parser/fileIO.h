/* 
 * File:   fileIO.h
 * Author: mmpi
 *
 * Created on September 26, 2013, 4:14 PM
 */

#ifndef FILEIO_H
#define	FILEIO_H

extern "C" {
  #include <sys/stat.h>
  #include <dirent.h>
  #include <errno.h>
}

#include <string>
#include <sstream>

inline bool fileExists(const std::string& fileName) {
  struct stat status;   
  return (stat(fileName.c_str(), &status) == 0); 
}

inline bool folderExists(const std::string& folderName) {
  struct stat status;   
  return (stat(folderName.c_str(), &status) == 0) && (status.st_mode & S_IFDIR); 
}

inline bool removeFolderContent(const std::string& folderName) {
  bool ret = true;
  struct dirent *nextFile;
  DIR *folder = opendir(folderName.c_str());
  while(nextFile = readdir(folder)) {
    // build the full path for each file in the folder
    std::string fileName(nextFile->d_name);
    if((fileName.compare(".")!=0) && (fileName.compare("..")!=0)) {
      std::stringstream ss;
      ss << folderName << "/" << fileName;
      if(0!=remove(ss.str().c_str())) {
        ret = false;
      }
    }
  }
  return ret;
}

inline bool recursivelyCreateFolder(std::string folderName) {
  size_t pre=0, pos;
  std::string dir;
  int mdret;

  if(folderName[folderName.size()-1]!='/'){
    folderName += '/';
  }

  while((pos=folderName.find_first_of('/', pre))!=std::string::npos) {
    dir = folderName.substr(0, pos++);
    pre = pos;
    if(dir.size()==0) continue; // if leading / first time is 0 length
    if((mdret=mkdir(dir.c_str(), 0755)) && errno!=EEXIST){
      return false;
    }
  }
  return true;
}


#endif	/* FILEIO_H */

