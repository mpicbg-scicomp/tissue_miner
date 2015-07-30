#include "Matrix2x2.h"

#include "Triangle.h"
#include "LogTriangleDeformation.h"

bool LogTriangleDeformation::startFile(const std::string &Path) {
  if(!LogFile::startFile(Path)) {
    return false;
  }
  // write header line
  _outputFile << "# ";
  _outputFile << "time (mid interval)" << HeaderSeparator;
  _outputFile << "delta t" << HeaderSeparator;
  logQuantityNames();
  _outputFile << std::endl;
  return true;
}


void LogTriangleDeformation::addDeformation(const TissueState *s, const TissueState *next) {
  prepareQuantities();
  for(unsigned int j=0; j<s->numberOfVertices(); ++j) {
    Vertex *v = s->vertex(j);
    for(unsigned int k=0; k<v->triangles.size(); ++k) {
      addTriangle(v->triangles[k]);
    }
  }
  _outputFile << 0.5*(s->time()+next->time()) << DataSeparator;
  _outputFile << next->time()-s->time() << DataSeparator;
  logQuantities();
  _outputFile << std::endl;
}


void LogTriangleDeformation::addTriangle(const Triangle *t) {
  const double AreaCutoff = 1e-10;
  
  if(!t->deformedTriangle) {
    std::cout << "TriangleDeformationAverage::addTriangle: no deformed triangle there!" << std::endl;
    throw std::exception();
  }
  // don't count in triangle, if area of deformed triangle below numerical cutoff
  if(t->deformedTriangle->area<AreaCutoff) {
    return;
  }
  // deformation matrix
  Matrix2x2 tU(t->deformedTriangle->t*t->t.inverse());
  // don't count in, if I can not decompose deformation matrix
  double deltaSU, deltaThetaU;
  Nematic2D deltaQU;
  if(!tU.decomposeIntoScalingRotationShear(deltaSU, deltaThetaU, deltaQU)) {
    return;
  }
  
  // OK, this is a triangle which I count in
  _sumArea1 += t->area;
  _wSumExp2DeltaSU += t->area * exp(2*deltaSU);
  _wSumDeltaThetaU += t->area * deltaThetaU;
  _wSumTUTilde += t->area * tU.symmetricTracelessPart();
  _wSumQAfterDeformation += t->area * t->deformedTriangle->Q;
  _wSumDeltaQU += t->area * deltaQU;
  _wSumQExp2DeltaSU += t->area * exp(2*deltaSU) * t->deformedTriangle->Q;
  _wSumQDeltaThetaU += t->area * deltaThetaU * 0.5*(t->Q + t->deformedTriangle->Q);
  double deltaPsiTheta = t->deformedTriangle->Q.angle()-t->Q.angle() -(t->deformedTriangle->theta-t->theta);
  deltaPsiTheta -= 2.0*M_PI*round(deltaPsiTheta/2.0/M_PI);
  Nematic2D Qbar(0.5*(t->Q + t->deformedTriangle->Q));
  double Q = Qbar.norm();
  double sinh2QOverQ = (Q==0)?2:sinh(2*Q)/Q;
  _wSumNonLinCorrection += t->area * deltaPsiTheta*(sinh2QOverQ-2*cosh(2*Q)) * Qbar;
}
void LogTriangleDeformation::prepareQuantities() {
  _sumArea1 = _wSumDeltaThetaU = _wSumExp2DeltaSU = 0;
  _wSumTUTilde.set(0,0);
  _wSumQAfterDeformation = _wSumDeltaQU = _wSumQExp2DeltaSU = _wSumQDeltaThetaU = _wSumNonLinCorrection = _wSumTUTilde;
}
void LogTriangleDeformation::logQuantities() {
  _outputFile << _sumArea1;
  _outputFile << DataSeparator << _wSumExp2DeltaSU/_sumArea1;
  _outputFile << DataSeparator << _wSumDeltaThetaU/_sumArea1;
  _wSumTUTilde /= _sumArea1;
  _outputFile << DataSeparator << _wSumTUTilde.c1() << DataSeparator << _wSumTUTilde.c2();
  _wSumQAfterDeformation /= _sumArea1;
  _outputFile << DataSeparator << _wSumQAfterDeformation.c1() << DataSeparator << _wSumQAfterDeformation.c2();
  _wSumDeltaQU /= _sumArea1;
  _outputFile << DataSeparator << _wSumDeltaQU.c1() << DataSeparator << _wSumDeltaQU.c2();
  _wSumQExp2DeltaSU /= _sumArea1;
  _outputFile << DataSeparator << _wSumQExp2DeltaSU.c1() << DataSeparator << _wSumQExp2DeltaSU.c2();
  _wSumQDeltaThetaU /= _sumArea1;
  _outputFile << DataSeparator << _wSumQDeltaThetaU.c1() << DataSeparator << _wSumQDeltaThetaU.c2();
  _wSumNonLinCorrection /= _sumArea1;
  _outputFile << DataSeparator << _wSumNonLinCorrection.c1() << DataSeparator << _wSumNonLinCorrection.c2();
}
void LogTriangleDeformation::logQuantityNames() {
  _outputFile << "sum area of the triangles counted for the averages" << HeaderSeparator;
  _outputFile << "<exp(2*Delta s_u)>_t" << HeaderSeparator;
  _outputFile << "<Delta theta_u>_t" << HeaderSeparator;
  _outputFile << "<tilde T_{u,1}>_t" << HeaderSeparator << "<tilde T_{u,2}>_t" << HeaderSeparator;
  _outputFile << "<Q_1 (after deformation)>_t" << HeaderSeparator << "<Q_2 (after deformation)>_t" << HeaderSeparator;
  _outputFile << "<Delta Q_{u,1}>_t" << HeaderSeparator << "<Delta Q_{u,2}>_t" << HeaderSeparator;
  _outputFile << "<Q_1*exp(2*Delta s_u)>_t" << HeaderSeparator << "<Q_2*exp(2*Delta s_u)>_t" << HeaderSeparator;
  _outputFile << "<bar Q_1*Delta theta_u>_t" << HeaderSeparator << "<bar Q_2*Delta theta_u>_t" << HeaderSeparator;
  _outputFile << "<(Delta psi-Delta theta)*(sinh(2Q)/Q-2cosh(q))*bar Q_1>_t" << HeaderSeparator << "<(Delta psi-Delta theta)*(sinh(2Q)/Q-2cosh(q))*bar Q_2>_t" << HeaderSeparator;
}
