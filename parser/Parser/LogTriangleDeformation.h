/* 
 * File:   TriangleDeformationAverage.h
 * Author: mmpi
 *
 * Created on 30. September 2013, 18:55
 */

#ifndef LOGTRIANGLEDEFORMATION_H
#define	LOGTRIANGLEDEFORMATION_H

#include "TissueState.h"
#include "Triangle.h"
#include "LogFile.h"

class LogTriangleDeformation : public LogFile {
public:
  virtual bool startFile(const std::string &Path);
  void addDeformation(const TissueState *s, const TissueState *next);

private:
  void logQuantityNames();
  void prepareQuantities();
  void addTriangle(const Triangle *t);
  void logQuantities();

  double _sumArea1, _wSumExp2DeltaSU, _wSumDeltaThetaU;
  Nematic2D _wSumTUTilde, _wSumQAfterDeformation, _wSumDeltaQU, _wSumQExp2DeltaSU, _wSumQDeltaThetaU, _wSumNonLinCorrection;
};

#endif	/* TRIANGLEDEFORMATIONAVERAGE_H */

