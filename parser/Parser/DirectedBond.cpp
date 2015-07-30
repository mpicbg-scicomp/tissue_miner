#include <algorithm>
#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif
#include "Triangle.h"
#include "TissueState.h"
#include "DirectedBond.h"

void DirectedBond::createTriangle() {
  if(bondTriangle) {
    delete bondTriangle;
  }
  if(deformedBond) {
    deformedBond->bondTriangle = new Triangle(deformedBond->rightVertex->r, deformedBond->leftVertex->r, deformedBond->cell->r);
    bondTriangle = new Triangle(rightVertex->r, leftVertex->r, cell->r, deformedBond->bondTriangle);
  } else {
    bondTriangle = new Triangle(rightVertex->r, leftVertex->r, cell->r);
  }
}

DirectedBond::~DirectedBond() {
  if(bondTriangle) {
    delete bondTriangle;
  }
}


#ifdef USE_NETCDF

void DirectedBond::DataFileIO::allocateBuffers(int numberOfFields) {
  _numberOfFields = numberOfFields;
  if(_numberOfFields>0) {
    _conjBond = new int[_numberOfFields];
    _cell = new int[_numberOfFields];
    _rightVertex = new int[_numberOfFields];
    _leftVertex = new int[_numberOfFields];
    _leftBondSeenFromVertex = new int[_numberOfFields];
    _leftBondSeenFromCell = new int[_numberOfFields];
  }
}
 
DirectedBond::DataFileIO::~DataFileIO() {
  if(_numberOfFields>0) {
    delete[] _conjBond;
    delete[] _cell;
    delete[] _rightVertex;
    delete[] _leftVertex;
    delete[] _leftBondSeenFromVertex;
    delete[] _leftBondSeenFromCell;
  }  
}
   
bool DirectedBond::DataFileIO::load(TissueState &s, const NcFile &dataFile) {
 // Create the cell dimension
  NcDim *bondsDim;
  if(!(bondsDim=dataFile.get_dim("bond index"))) {
    return true;
  }
  allocateBuffers(bondsDim->size());

  // Create the variables
  if(!(_conjBondVar = dataFile.get_var("DirectedBond::conjBond"))) return false;
  if(!(_cellVar = dataFile.get_var("DirectedBond::cell"))) return false;
  if(!(_rightVertexVar = dataFile.get_var("DirectedBond::rightVertex"))) return false;
  if(!(_leftVertexVar = dataFile.get_var("DirectedBond::leftVertex"))) return false;
  if(!(_leftBondSeenFromVertexVar = dataFile.get_var("DirectedBond::_leftBondSeenFromVertex"))) return false;
  if(!(_leftBondSeenFromCellVar = dataFile.get_var("DirectedBond::_leftBondSeenFromCell"))) return false;
 
  // Create arrays and load the data into them
  if(!_conjBondVar->get(&_conjBond[0], _numberOfFields)) return false;
  if(!_cellVar->get(&_cell[0], _numberOfFields)) return false;
  if(!_rightVertexVar->get(&_rightVertex[0], _numberOfFields)) return false;
  if(!_leftVertexVar->get(&_leftVertex[0], _numberOfFields)) return false;
  if(!_leftBondSeenFromVertexVar->get(&_leftBondSeenFromVertex[0], _numberOfFields)) return false;
  if(!_leftBondSeenFromCellVar->get(&_leftBondSeenFromCell[0], _numberOfFields)) return false;

  // keep map tid -> bond pointer
  std::vector<DirectedBond*> tidToBond;
  // Create the bonds
  for(int i=0; i<_numberOfFields; ++i) {
    DirectedBond *b = s.newBond();
    tidToBond.push_back(b);
    b->cell = s.cell(_cell[i]);
    b->cell->bonds.push_back(b);
    b->rightVertex = s.vertex(_rightVertex[i]);
    b->rightVertex->bonds.push_back(b);
    b->leftVertex = s.vertex(_leftVertex[i]);
  }
  // connect conjugated bonds
  std::map<DirectedBond*,DirectedBond*> leftBondSeenFromVertexMap, leftBondSeenFromCellMap;
  for(int i=0; i<_numberOfFields; ++i) {
    DirectedBond *b = tidToBond[i];
    if(_conjBond[i]<0) {
      b->conjBond = NULL;
    } else {
      b->conjBond = tidToBond[_conjBond[i]];
    }
    leftBondSeenFromVertexMap[b] = tidToBond[_leftBondSeenFromVertex[i]];
    leftBondSeenFromCellMap[b] = tidToBond[_leftBondSeenFromCell[i]];
  }
  // sort bonds around vertices
  for(int i=0; i<s.numberOfVertices(); ++i) {
    Vertex *v = s.vertex(i);
    for(unsigned int j=1; j<v->bonds.size()-1; ++j) {
      DirectedBond *next = leftBondSeenFromVertexMap[v->bonds[j-1]];
      bool found = false;
      for(unsigned int k=j; k<v->bonds.size(); ++k) {
        if(v->bonds[k]==next) {
          if(k!=j) {
            v->bonds[k] = v->bonds[j];
            v->bonds[j] = next;
          }
          found = true;
          break;
        }
      }
      if(!found) {
        std::cout << "DirectedBond::DataFileIO::load: Did not find next bond for vertex!" << std::endl;
        return false;
      }
    }
  }  
  // sort bonds around cells
  for(TissueState::CellIterator it=s.beginCellIterator(); it!=s.endCellIterator(); ++it) {
    Cell *c = s.cell(it);
    for(unsigned int j=1; j<c->bonds.size()-1; ++j) {
      DirectedBond *next = leftBondSeenFromCellMap[c->bonds[j-1]];
      bool found = false;
      for(unsigned int k=j; k<c->bonds.size(); ++k) {
        if(c->bonds[k]==next) {
          if(k!=j) {
            c->bonds[k] = c->bonds[j];
            c->bonds[j] = next;
          }
          found = true;
          break;
        }
      }
      if(!found) {
        std::cout << "DirectedBond::DataFileIO::load: Did not find next bond for cell!" << std::endl;
        return false;
      }
    }
  }  
  
  return true;    
}
  
bool DirectedBond::DataFileIO::save(NcFile &dataFile, const TissueState &s) {
  // This condition is to allow for empty networks to be saved.
  // Some netCDF versions have problems with dimension length 0. They think, it's the unlimited dimension.
  // So, if there are no vertices, nothing will be written to the file. Not even the dimension.
  if(s.numberOfBonds()) {
    // Create the bond dimension
    NcDim *bondsDim;
    if(!(bondsDim = dataFile.add_dim("bond index", s.numberOfBonds()))) return false;

    // Create the variables
    if(!(_conjBondVar = dataFile.add_var("DirectedBond::conjBond", ncInt, bondsDim))) return false;
    if(!(_cellVar = dataFile.add_var("DirectedBond::cell", ncInt, bondsDim))) return false;
    if(!(_rightVertexVar = dataFile.add_var("DirectedBond::rightVertex", ncInt, bondsDim))) return false;
    if(!(_leftVertexVar = dataFile.add_var("DirectedBond::leftVertex", ncInt, bondsDim))) return false;
    if(!(_leftBondSeenFromVertexVar = dataFile.add_var("DirectedBond::_leftBondSeenFromVertex", ncInt, bondsDim))) return false;
    if(!(_leftBondSeenFromCellVar = dataFile.add_var("DirectedBond::_leftBondSeenFromCell", ncInt, bondsDim))) return false;
 
    // generate map: vertex pointer -> tid
    std::map<Vertex*,int> vertexPointerToTid;
    std::map<DirectedBond*,DirectedBond*> leftBondSeenFromVertexMap;
    for(int i=0; i<s.numberOfVertices(); ++i) {
      Vertex *v = s.vertex(i);
      vertexPointerToTid[v] = i;
      for(unsigned int j=0; j<v->bonds.size(); ++j) {
        leftBondSeenFromVertexMap[v->bonds[j]] = v->bonds[(j+1)%v->bonds.size()];
      }
    }
    // generate map: bond pointer -> tid
    std::map<DirectedBond*,int> bondPointerToTid;
    int i = 0;
    for(TissueState::BondConstIterator it=s.beginBondIterator(); it!=s.endBondIterator(); ++it) {
      bondPointerToTid[s.bond(it)] = i;
      ++i;
    }
    // generate map for bond sorting within cells
    std::map<DirectedBond*,DirectedBond*> leftBondSeenFromCellMap;
    for(TissueState::CellConstIterator it=s.beginCellIterator(); it!=s.endCellIterator(); ++it) {
      Cell *c = s.cell(it);
      for(unsigned int j=0; j<c->bonds.size(); ++j) {
        leftBondSeenFromCellMap[c->bonds[j]] = c->bonds[(j+1)%c->bonds.size()];
      }
    }

    // fill the buffers
    allocateBuffers(s.numberOfBonds());
    i = 0;
    for(TissueState::BondConstIterator it=s.beginBondIterator(); it!=s.endBondIterator(); ++it) {
      DirectedBond *b = s.bond(it);
      if(b->conjBond) {
        _conjBond[i] = bondPointerToTid.at(b->conjBond);
      } else {
        _conjBond[i] = -1;
      }
      _cell[i] = b->cell->id;
      _rightVertex[i] = vertexPointerToTid.at(b->rightVertex);
      _leftVertex[i] = vertexPointerToTid.at(b->leftVertex);
      _leftBondSeenFromVertex[i] = bondPointerToTid.at(leftBondSeenFromVertexMap.at(b));
      _leftBondSeenFromCell[i] = bondPointerToTid.at(leftBondSeenFromCellMap.at(b));
      ++i;
    }
    
    // put the records
    if(!_conjBondVar->put(_conjBond, s.numberOfBonds())) return false;
    if(!_cellVar->put(_cell, s.numberOfBonds())) return false;
    if(!_rightVertexVar->put(_rightVertex, s.numberOfBonds())) return false;
    if(!_leftVertexVar->put(_leftVertex, s.numberOfBonds())) return false;
    if(!_leftBondSeenFromVertexVar->put(_leftBondSeenFromVertex, s.numberOfBonds())) return false;
    if(!_leftBondSeenFromCellVar->put(_leftBondSeenFromCell, s.numberOfBonds())) return false;
  }
  return true;  
}

#endif  
