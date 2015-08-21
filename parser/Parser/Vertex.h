/* 
 * File:   Vertex.h
 * Author: mmpi
 *
 * Created on September 13, 2013, 10:58 AM
 */

#ifndef VERTEX_H
#define	VERTEX_H

#include <vector>

#include "Object2D.h"
#include "DirectedBond.h"

class TissueState;

/** Contains Vertex properties and references to bonds and cells */
struct Vertex {
  Vertex(const Vector2D &v) : r(v) {}
  Vertex(const Vertex &other) : r(other.r) {}
  
  /** These are pointers to all bonds of this vertex sorted in ccw order. The cell a directed bond points to is right of the bonds seen from this vertex. */
  std::vector<DirectedBond *> bonds;
  
  /** position of the vertex */
  Vector2D r;
};

#endif	/* VERTEX_H */

