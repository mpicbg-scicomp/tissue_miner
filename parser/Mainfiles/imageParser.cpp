/* 
 * File:   wingDeformation.cpp
 * Author: mmpi
 *
 * Created on September 13, 2013, 10:48 AM
 */

#include <stdlib.h>
#include <iostream>

#include "MovieData.h"
#include "DeformingTissue.h"


int main(int argc, char** argv) {
  if(argc<4) {
    cout << PROGRAMNAME << ", Version: " << VERSION << endl;
    cout << "Usage: ./imageParser <root path> <movie folder name> <frame number mask>" << endl;
    cout << endl;
    cout << "Example: ./imageParser /movieSegmentation WT_111102 %03d" << endl << endl;
    return 0;
  }
  const MovieData Movie(argv[1], argv[2], argv[3]);
  
  DeformingTissue t(Movie);
  if(!t.loadFromSegmentedData()) {
    return 1;
  }         
  if(!t.exportToDbTables()) {
    return 1;
  }

  
  // check number of daughters with same ids
  multiset<CellIndex> daughters;
  for(int i=0; i<t.numberOfFrames(); ++i) {
    TissueState *s=t.frame(i);
    for(TissueState::CellConstIterator it=s->beginCellIterator(); it!=s->endCellIterator(); ++it) {
      Cell *c = s->cell(it);
//      if((c->duringTransitionBefore==Cell::Divided) && (c->mother!=c->id)) {
//        daughters.insert(c->id);
//      }
      if(c->daughter) {
        daughters.insert(c->daughter);
      }
    }
  }
  int count(0);
  CellIndex example(0);
  for(multiset<CellIndex>::const_iterator it=daughters.begin(); it!=daughters.end(); ++it) {
    if(daughters.count(*it)>1) {
      if(!count) {
        example = *it;
      }
      ++count;
    }
  }
  cout << count << " cell ids appeared more than once as daughter cell ids." << endl;
  if(count) {
    cout << "For example daughter cell id: " << example << endl;
  }
    
  cout << "Done." << endl;

  return 0;
}

