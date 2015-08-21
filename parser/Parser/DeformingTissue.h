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
#include "TissueState.h"


class DeformingTissue {
public:
  /** creates empty frame set */
  DeformingTissue(const MovieData &movie) : _movie(movie) {};
  ~DeformingTissue();

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

  /** exports the whole topology data etc. to text files for database */
  bool exportToDbTables() const;
  
  /** number of frames */
  int numberOfFrames() const { return _frames.size(); }
  /** tissue state */
  TissueState *frame(const int index) const { return _frames[index]; }

private:
  const MovieData &_movie;
  const static std::string _dbSubFolder;
  
  std::vector<TissueState*> _frames;
  
  /** loads the data from the divisions file belonging to frame i */
  bool parseDivisionsFile(const string TrackedCellsFileName, const string DivisionsFileName, const int frame);
  /** db folder name, from movie and label */
  std::string dbFolderName() const;
  /** creates empty state, but keeps connection to movie */
  void cleanUp();
};

#endif	/* DEFORMINGTISSUE_H */

