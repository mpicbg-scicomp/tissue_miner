#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include "Color.h"
#include "Object2D.h"
#include "MovieData.h"

using namespace std;

const string MovieData::MovieInformationFileName("movieInformation.dat");
const string MovieData::TimesInSecondsFileName("cumultimesec.txt");
const string MovieData::SegmentedFrameSubFolderMask("Segmentation");
const string MovieData::OriginalImageFileName("original.png");
const string MovieData::TrackedCellsImageFileName("tracked_cells_resized.tif");
const string MovieData::DaughterCellsImageFileName("dividing_cells.tif");
const string MovieData::ApoptosisImageFileName("apoptotic_cells.png");

 
MovieData::MovieData(const string &root_, const string &movieName_, const string &frameMask_) :  root(root_), movieName(movieName_), frameMask(frameMask_), moviePath(generateMoviePath()) {
  stringstream miFileName("");

  // get movie information
  miFileName << moviePath << SegmentedFrameSubFolderMask << "/" << MovieInformationFileName;
  ifstream movieInformationFile(miFileName.str().c_str());
  if(movieInformationFile.fail()) {
    cout << "MovieData::MovieData: Unable to open " << miFileName.str() << "!" << endl;
    throw exception();
  }
  int cti;
  movieInformationFile >> width >> height >> numFrames >> firstTimeInHApf >> timeIntervalInH >> cti >> resolutionInUmPerPixel;
  if(movieInformationFile.fail()) {
    cout << "MovieData::MovieData: " << miFileName.str() << "has invalid format!" << endl;
    throw exception();
  }
  movieInformationFile.close();
  constTimeInterval = (cti!=0);
  
  // set times (in hAPF):
  frameTime = new double[numFrames];
  if(constTimeInterval) {
    for(int frame=0; frame<numFrames; ++frame) {
      frameTime[frame] = firstTimeInHApf+timeIntervalInH*frame;
    }
  } else {
    stringstream tisFileName("");
    tisFileName << moviePath << SegmentedFrameSubFolderMask << "/" << TimesInSecondsFileName;
    ifstream timesInSecondsFile(tisFileName.str().c_str());
    if(timesInSecondsFile.fail()) {
      cout << "MovieData::MovieData: Unable to open " << tisFileName.str() << "!" << endl;
      throw exception();
    }
    for(int frame=0; frame<numFrames; ++frame) {
      double tis;
      timesInSecondsFile >> tis;
      if(timesInSecondsFile.fail()) {
        cout << "MovieData::MovieData: " << tisFileName.str() << "has invalid format or too few times!" << endl;
        throw exception();
      }
      frameTime[frame] = firstTimeInHApf+timeIntervalInH*tis;
    }
    timesInSecondsFile.close();
  }
}

MovieData::~MovieData() {
  delete frameTime;
}

string MovieData::generateMoviePath() const {
  stringstream folder("");
  folder << root << "/" << movieName << "/";
  return folder.str();
}

string MovieData::segmentedFrameFolder(const int Frame) const {
  char frameSubFolder[100];
  sprintf(frameSubFolder, frameMask.c_str(), Frame);
  stringstream folder("");
  folder << moviePath << "/" << SegmentedFrameSubFolderMask << "/" << movieName << "_" << frameSubFolder << "/";
  return folder.str();
}

string MovieData::movieFilePath(const string &FileName) const {
  stringstream fileName("");
  fileName << moviePath << FileName;
  return fileName.str();
}

string MovieData::segmentedFrameFilePath(const int Frame, const string &FileName) const {
  stringstream fileName("");
  fileName << segmentedFrameFolder(Frame) << FileName;
  return fileName.str();
}

string MovieData::originalImageFilePath(const int Frame) const {
  return segmentedFrameFilePath(Frame, OriginalImageFileName);
}

string MovieData::trackedCellsImageFilePath(const int Frame) const {
  return segmentedFrameFilePath(Frame, TrackedCellsImageFileName);
}

string MovieData::daughterCellsImageFilePath(const int Frame) const {
  return segmentedFrameFilePath(Frame, DaughterCellsImageFileName);
}

