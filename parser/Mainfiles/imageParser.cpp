/* 
 * File:   wingDeformation.cpp
 * Author: mmpi
 *
 * Created on September 13, 2013, 10:48 AM
 */

#include <stdlib.h>

#include "MovieData.h"
#include "QtRasterImage.h"
#include "QtVectorImage.h"
#include "DeformingTissue.h"

void drawNetwork(QtImage &img, const IsotropicTransformation &Ref2Pixel, const TissueState &s) {
  // add offset because of the way, pdfs are stored
  s.drawBonds(img, IsotropicTransformation(Vector2D(0.5, -0.5))*Ref2Pixel, Black, 0.5, true);
  s.drawVertices(img, IsotropicTransformation(Vector2D(0.5, -0.5))*Ref2Pixel, Red, 0.5);
}

void drawDivisionPairs(QtImage &img, const IsotropicTransformation &Ref2Pixel, const TissueState &s) {
  s.drawTime(img);
  s.drawDivisionPairs(img, Ref2Pixel, Red, Green, 8, Blue, 3);
}

void drawTransitionBefore(QtImage &img, const IsotropicTransformation &Ref2Pixel, const TissueState &s) {
  s.drawTime(img);
  s.drawCellsTransitionBefore(img, Ref2Pixel, 8);
}

void drawTransitionAfter(QtImage &img, const IsotropicTransformation &Ref2Pixel, const TissueState &s) {
  s.drawTime(img);
  s.drawCellsTransitionAfter(img, Ref2Pixel, 8);
}

void drawPolarity(QtImage &img, const IsotropicTransformation &Ref2Pixel, const TissueState &s) {
  s.drawTime(img);
  s.drawCellsPolarityNormalizedByIntIntensityR(img, Ref2Pixel, Green, 100);
}

void drawPolarityBoxesNematicOrder(QtImage &img, const IsotropicTransformation &Ref2Pixel, const TissueState &s) {
  const double BoxSize = 300, LenFactor=1.0, Width=0.07*BoxSize;
  const int Nx = ceil(img.width()/BoxSize), Ny = ceil(img.height()/BoxSize);
  Nematic2D sumPolarities[Nx][Ny];
  double sumNorms[Nx][Ny];
  int count[Nx][Ny];
  
  for(int i=0; i<Nx; ++i) {
    for(int j=0; j<Ny; ++j) {
      sumPolarities[i][j].set(0,0);
      sumNorms[i][j] = 0;
      count[i][j] = 0;
    }
  }
  for(TissueState::CellConstIterator it=s.beginCellIterator(); it!=s.endCellIterator(); ++it) {
    Cell *c=s.cell(it);
    Vector2D r(Ref2Pixel.map(c->r));
    int i = int(r.x()/BoxSize);
    int j = int(r.y()/BoxSize);
    sumPolarities[i][j] += c->polarityR;
    sumNorms[i][j] += c->polarityR.norm();
    count[i][j] += 1;
  }
  s.drawTime(img);
  for(int i=0; i<Nx; ++i) {
    for(int j=0; j<Ny; ++j) {
      if(count[i][j]) {
        img.drawBar(Vector2D(BoxSize*(i+0.5), BoxSize*(j+0.5)), rotatedObject(LenFactor*BoxSize*sumPolarities[i][j]/sumNorms[i][j], Ref2Pixel.angle()), Green, Width);
      }
    }
  }
}

void drawPolarityBoxesNematicOrderIntensityWeight(QtImage &img, const IsotropicTransformation &Ref2Pixel, const TissueState &s) {
  const double BoxSize = 300, LenFactor=1.0, Width=0.07*BoxSize;
  const int Nx = ceil(img.width()/BoxSize), Ny = ceil(img.height()/BoxSize);
  Nematic2D sumPolarities[Nx][Ny];
  double sumNorms[Nx][Ny];
  int count[Nx][Ny];
  
  for(int i=0; i<Nx; ++i) {
    for(int j=0; j<Ny; ++j) {
      sumPolarities[i][j].set(0,0);
      sumNorms[i][j] = 0;
      count[i][j] = 0;
    }
  }
  for(TissueState::CellConstIterator it=s.beginCellIterator(); it!=s.endCellIterator(); ++it) {
    Cell *c=s.cell(it);
    Vector2D r(Ref2Pixel.map(c->r));
    int i = int(r.x()/BoxSize);
    int j = int(r.y()/BoxSize);
    sumPolarities[i][j] += c->polarityR/c->intIntensityR;
    sumNorms[i][j] += c->polarityR.norm()/c->intIntensityR;
    count[i][j] += 1;
  }
  s.drawTime(img);
  for(int i=0; i<Nx; ++i) {
    for(int j=0; j<Ny; ++j) {
      if(count[i][j]) {
        img.drawBar(Vector2D(BoxSize*(i+0.5), BoxSize*(j+0.5)), rotatedObject(LenFactor*BoxSize*sumPolarities[i][j]/sumNorms[i][j], Ref2Pixel.angle()), Green, Width);
      }
    }
  }
}

void drawPolarityBoxes(QtImage &img, const IsotropicTransformation &Ref2Pixel, const TissueState &s) {
  const double BoxSize = 300, LenFactor=20.0, Width=0.07*BoxSize;;
  const int Nx = ceil(img.width()/BoxSize), Ny = ceil(img.height()/BoxSize);
  Nematic2D sumPolarities[Nx][Ny];
  double sumIntensities[Nx][Ny];
  int count[Nx][Ny];
  
  for(int i=0; i<Nx; ++i) {
    for(int j=0; j<Ny; ++j) {
      sumPolarities[i][j].set(0,0);
      sumIntensities[i][j] = 0;
      count[i][j] = 0;
    }
  }
  for(TissueState::CellConstIterator it=s.beginCellIterator(); it!=s.endCellIterator(); ++it) {
    Cell *c=s.cell(it);
    Vector2D r(Ref2Pixel.map(c->r));
    int i = int(r.x()/BoxSize);
    int j = int(r.y()/BoxSize);
    sumPolarities[i][j] += c->polarityR;
    sumIntensities[i][j] += c->intIntensityR;
    count[i][j] += 1;
  }
  s.drawTime(img);
  for(int i=0; i<Nx; ++i) {
    for(int j=0; j<Ny; ++j) {
      if(count[i][j]) {
        img.drawBar(Vector2D(BoxSize*(i+0.5), BoxSize*(j+0.5)), rotatedObject(LenFactor*BoxSize*sumPolarities[i][j]/sumIntensities[i][j], Ref2Pixel.angle()), Green, Width);
      }
    }
  }
}

int main(int argc, char** argv) {
  QApplication app(argc, argv);
  
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
//  if(!t.save()) {
//    return 1;
//  }
//  if(!t.load("rawIncludingCellDivisions")) {
//    return 1;
//  }
  if(!t.exportToDbTables()) {
    return 1;
  }
//  t.removeCellsMiom();
//  const string MaskTag("blade");
//  const string MaskTag("hinge");
//  if(!t.applyMaskToFinalState(MaskTag)) {
//    return -1;
//  }
//  if(!t.save()) {
//    return 1;
//  }
//  if(!t.load(MaskTag)) {
//    return 1;
//  }
//  if(!t.writeDeformationData()) {
//    return 1;
//  }  
//  if(!t.exportAllCellIds()) {
//    return 1;
//  }  
//  if(!t.drawFrames("divisionPairs", drawDivisionPairs)) {
//    return 1;
//  }
//  if(!t.drawFrames("transitionBefore", drawTransitionBefore)) {
//    return 1;
//  }
//  if(!t.drawFrames("transitionAfter", drawTransitionAfter)) {
//    return 1;
//  }
//  if(!t.drawFrames("polarityBoxes_nematicOrder(300)", drawPolarityBoxesNematicOrder)) {
//    return 1;
//  }
//  if(!t.drawFrames("polarityBoxes_intIntensityNormalization(300)", drawPolarityBoxes)) {
//    return 1;
//  }
//  if(!t.drawFrames("polarityBoxes_nematicOrderIntensityWeight(300)", drawPolarityBoxesNematicOrderIntensityWeight)) {
//    return 1;
//  }
//  if(!t.drawFrames("polarity", drawPolarity, false)) {
//    return 1;
//  }
  
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

