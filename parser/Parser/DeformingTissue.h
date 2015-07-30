/* 
 * File:   DeformingTissue.h
 * Author: mmpi
 *
 * Created on September 25, 2013, 4:48 PM
 */

#ifndef DEFORMINGTISSUE_H
#define	DEFORMINGTISSUE_H

#include <vector>
#include <string>
#include "MovieData.h"
#include "QtImage.h"
#include "TissueState.h"


class DeformingTissue {
public:
  /** creates empty frame set */
  DeformingTissue(const MovieData &movie) : _movie(movie), _label("void") {};
  ~DeformingTissue();

#ifdef USE_NETCDF
  /** creates empty frame set */
  DeformingTissue(const MovieData &movie, const std::string label, const int maxFrames=-1);
  /** loads frame set data */
  bool load(const std::string label, const int maxFrames=-1);
  /** save frame set data */
  bool save() const;
#endif  

  /**
   * loads frames from tracked cells files from segmented movie
   * @param Movie
   * @param TrackedCellsImageMask
   * @return true, if successful
   */
  bool loadFromSegmentedData(const int MaxFrames=-1);
  
  /**
   * loads cell division data from whitespace-separated text file
   * @param FileName
   * @param skipRows
   * @return true if successful
   */
  bool loadCellDivisionData(const std::string FileName, const int skipRows=1);

  /** updates status of all cells */
  void fixCellStatus();

  /** removes cells moving into or out of the margin */
  void removeCellsMiom();
  
  /** selects region by applying mask image to last frame; only cells in non-black regions are kept; the red channel is used. */
  bool applyMaskToFinalState(const std::string MaskImageFileName);
  
  /** writes deformation data files */
  bool writeDeformationData();
  
  /** exports the whole topology data etc. to text files for database */
  bool exportToDbTables() const;
  
  /** exports all cell ids */
  bool exportAllCellIds() const;
  
  /** number of frames */
  int numberOfFrames() const { return _frames.size(); }
  /** tissue state */
  TissueState *frame(const int index) const { return _frames[index]; }

  /** draws series of tissue state images with original images as background */
  bool drawFrames(const std::string Tag, void (*drawFrame)(QtImage &img, const IsotropicTransformation &Ref2Pixel, const TissueState& s), const bool smaller=true) const;

  /** draws series of tissue state vector images with tracked cells images as background */
  bool drawFramesTcVec(const std::string Tag, void (*drawFrame)(QtImage &img, const IsotropicTransformation &Ref2Pixel, const TissueState& s)) const;
  
  /** draws the series of tissue states before division on vector images with tracked cells images as background */
  bool drawStatesBeforeDivisionDebug() const;
  
private:
  const MovieData &_movie;
  const static std::string _dbSubFolder;
  const static std::string _dataSubFolder;
  const static std::string _deformationFolder;
  const static std::string _imagesSubFolder;
  
  std::string _label;
  std::vector<TissueState*> _frames, _statesBeforeDivision, _statesAfterDeformation;
  
  /** loads the data from the divisions file belonging to frame i */
  bool parseDivisionsFile(const string TrackedCellsFileName, const string DivisionsFileName, const int frame);
  /** removes cells from all frames including cells appearing/disappearing at the margin */
  void removeCellsMiomHelper(std::set<CellIndex> &toBeRemoved);
  /** deformation data folder name, from movie and label */
  std::string textDataFolderName() const;
  /** file name of deformation data, from movie and label */
  std::string textDataFileName(const std::string &FileName) const;
  /** image folder name, from movie and label */
  std::string imageFolderName(const std::string &Tag) const;
  /** image file name, from movie and label */
  std::string imageFileName(const std::string &Tag, const std::string &Extension, const int frame) const;
  /** data folder name, from movie and label */
  std::string dataFolderName() const;
  /** data file name, from movie and label */
  std::string dataFileName(const int frame) const;
  /** db folder name, from movie and label */
  std::string dbFolderName() const;
  /** creates empty state, but keeps connection to movie */
  void cleanUp();
};

#endif	/* DEFORMINGTISSUE_H */

