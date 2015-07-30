/* 
 * File:   Vertex.h
 * Author: mmpi
 *
 * Created on September 13, 2013, 10:58 AM
 */

#ifndef VERTEX_H
#define	VERTEX_H

#include <vector>

#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif

#include "Object2D.h"
#include "DirectedBond.h"
#include "Triangle.h"

class TissueState;

/** Contains Vertex properties and references to bonds and cells */
struct Vertex {
  Vertex(const Vector2D &v) : r(v), deformedVertex(NULL), margin(false) {}
  Vertex(const Vertex &other) : r(other.r), deformedVertex(NULL) {}
  void createTriangles();
  void cleanUpTriangles();
  ~Vertex();
  
#ifdef USE_NETCDF
  /** for saving the data contained in all vertex objects */
  class DataFileIO {
  public:
    DataFileIO() : _numberOfFields(0) {}
    ~DataFileIO();
    bool load(TissueState &s, const NcFile &dataFile);
    bool save(NcFile &dataFile, const TissueState &s);
  private:
    void allocateBuffers(int numberOfFields);
    int _numberOfFields;
    NcVar *_xVar, *_yVar;
    double *_x, *_y;
  };
#endif  
  
  /** These are pointers to all bonds of this vertex sorted in ccw order. The cell a directed bond points to is right of the bonds seen from this vertex. */
  std::vector<DirectedBond *> bonds;
  
  /** position of the vertex */
  Vector2D r;
  
  /** is a vertex directly at the margin */
  bool margin;
  
  /** deformed vertex */
  Vertex *deformedVertex;
  
  /** list of triangles associated with this vertex */
  std::vector<Triangle *> triangles;
};

#endif	/* VERTEX_H */

