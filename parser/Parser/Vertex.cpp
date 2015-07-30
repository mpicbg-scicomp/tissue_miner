#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif
#include "TissueState.h"
#include "Vertex.h"

void Vertex::createTriangles() {
  // NOTE: there are several alternatives to define the triangles. This is the only place which is affected by this choice
  // but: vertex positions need to be updated for intermediate frames!
  cleanUpTriangles();
  
  // check for margin vertex, first
  bool margin = false;
  for(unsigned int j=0; j<bonds.size(); ++j) {
    if(!bonds[j]->conjBond) {
      margin = true;
    }
  }
  
  if(!margin) {
    if(deformedVertex) {
      if(deformedVertex->triangles.size()>0) {
        std::cout << "Vertex::createTriangles: deformed vertex already has triangles defined!" << std::endl;
        throw std::exception();
      }
      if(deformedVertex->bonds.size()!=bonds.size()) {
        std::cout << "Vertex::createTriangles: deformed vertex has different bond number!" << std::endl;
        throw std::exception();
      }
      if(bonds.size()==3) {
        Triangle *t = new Triangle(deformedVertex->bonds[0]->cell->r, deformedVertex->bonds[1]->cell->r, deformedVertex->bonds[2]->cell->r);
        deformedVertex->triangles.push_back(t);
        triangles.push_back(new Triangle(bonds[0]->cell->r, bonds[1]->cell->r, bonds[2]->cell->r, t));
      } else if(bonds.size()>3) {
        // center of all surrounding cells
        Vector2D center(0,0), centerDeformed(0,0);
        for(unsigned int j=0; j<bonds.size(); ++j) {
          center += bonds[j]->cell->r;
          centerDeformed += bonds[j]->cell->deformedCell->r;
        }
        center /= bonds.size();
        centerDeformed /= bonds.size();

        // create triangles
        Cell *lastCell = bonds[bonds.size()-1]->cell;
        for(unsigned int j=0; j<bonds.size(); ++j) {
          Cell *c = bonds[j]->cell;
          Triangle *t = new Triangle(centerDeformed, lastCell->deformedCell->r, c->deformedCell->r);
          deformedVertex->triangles.push_back(t);
          triangles.push_back(new Triangle(center, lastCell->r, c->r, t));
          lastCell = c;
        }
      } else {
        std::cout << "Vertex::createTriangles: non margin vertex has less than 3 bonds!" << std::endl;
        throw std::exception();
      }
    } else {
      // !deformedVertex
      if(bonds.size()==3) {
        triangles.push_back(new Triangle(bonds[0]->cell->r, bonds[1]->cell->r, bonds[2]->cell->r));
      } else if(bonds.size()>3) {
        // center of all surrounding cells
        Vector2D center(0,0);
        for(unsigned int j=0; j<bonds.size(); ++j) {
          center += bonds[j]->cell->r;
        }
        center /= bonds.size();

        // create triangles
        Cell *lastCell = bonds[bonds.size()-1]->cell;
        for(unsigned int j=0; j<bonds.size(); ++j) {
          Cell *c = bonds[j]->cell;
          triangles.push_back(new Triangle(center, lastCell->r, c->r));
          lastCell = c;
        }
      } else {
        std::cout << "Vertex::createTriangles: non margin vertex has less than 3 bonds!" << std::endl;
        throw std::exception();
      }
    }
  } else {
    // do nothing
  }
}

void Vertex::cleanUpTriangles() {
  for(unsigned int i=0; i<triangles.size(); ++i) {
    delete triangles[i];
  }
}

Vertex::~Vertex() {
  cleanUpTriangles();
}

#ifdef USE_NETCDF

void Vertex::DataFileIO::allocateBuffers(int numberOfFields) {
  _numberOfFields = numberOfFields;
  if(_numberOfFields>0) {
    _x = new double[_numberOfFields];
    _y = new double[_numberOfFields];
  }
}

Vertex::DataFileIO::~DataFileIO() {
  if(_numberOfFields>0) {
    delete[] _x;
    delete[] _y;
  }  
}
   
bool Vertex::DataFileIO::load(TissueState &s, const NcFile &dataFile) {
 // Create the vertex dimension
  NcDim *verticesDim;
  if(!(verticesDim=dataFile.get_dim("vertex index"))) {
    return true;
  }
  allocateBuffers(verticesDim->size());
  
  // Create the variables
  if(!(_xVar = dataFile.get_var("Vertex::x"))) return false;
  if(!(_yVar = dataFile.get_var("Vertex::y"))) return false;

  // Create arrays and load the data into them
  if(!_xVar->get(&_x[0], _numberOfFields)) return false;
  if(!_yVar->get(&_y[0], _numberOfFields)) return false;

  // Create the vertices
  for(int i=0; i<_numberOfFields; ++i) {
    Vector2D r(_x[i], _y[i]);
    s.newVertex(r);
  }
  
  return true;  
}

bool Vertex::DataFileIO::save(NcFile &dataFile, const TissueState &s) {
  // This condition is to allow for empty networks to be saved.
  // Some netCDF versions have problems with dimension length 0. They think, it's the unlimited dimension.
  // So, if there are no vertices, nothing will be written to the file. Not even the dimension.
  if(s.numberOfVertices()) {
    // Create the bond dimension
    NcDim *verticesDim;
    if(!(verticesDim = dataFile.add_dim("vertex index", s.numberOfVertices()))) return false;

    // Create the variables
    if(!(_xVar = dataFile.add_var("Vertex::x", ncDouble, verticesDim))) return false;
    if(!(_yVar = dataFile.add_var("Vertex::y", ncDouble, verticesDim))) return false;
    
    // fill the buffers
    allocateBuffers(s.numberOfVertices());
    for(int i=0; i<s.numberOfVertices(); ++i) {
      Vector2D &v = s.vertex(i)->r;
      _x[i] = v.x();
      _y[i] = v.y();
    }
    
    // put the records
    if(!_xVar->put(_x, s.numberOfVertices())) return false;
    if(!_yVar->put(_y, s.numberOfVertices())) return false;
  }
  return true;  
}

#endif
