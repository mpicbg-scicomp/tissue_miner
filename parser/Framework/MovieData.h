#ifndef MOVIE_DATA_H
#define MOVIE_DATA_H

#include <limits>
#include <string>
#include "Color.h"
#include "Object2D.h"

using namespace std;

struct MovieData {
  MovieData(const string &root_, const string &movieName_, const string &frameMask_);
  ~MovieData();

  string movieFilePath(const string &FileName) const;
  string segmentedFrameFilePath(const int Frame, const string &FileName) const;
  string originalImageFilePath(const int Frame) const;
  string trackedCellsImageFilePath(const int Frame) const;
  string daughterCellsImageFilePath(const int Frame) const;
  
  const string root;
  const string movieName;
  const string frameMask;
  const string moviePath;
  
  // content of movie information file (in this order)
  int width, height, numFrames;
  double firstTimeInHApf;
  double timeIntervalInH;
  bool constTimeInterval;
  double resolutionInUmPerPixel;

  // time points:
  double *frameTime;
  
private:
  string generateMoviePath() const;
  string segmentedFrameFolder(const int Frame) const;

  const static string SegmentedFrameSubFolderMask, MovieInformationFileName, TimesInSecondsFileName, OriginalImageFileName, TrackedCellsImageFileName, DaughterCellsImageFileName, ApoptosisImageFileName;
};

#endif
