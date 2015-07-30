#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif
#include "TissueState.h"
#include "Cell.h"

#ifdef USE_NETCDF
void Cell::DataFileIO::allocateBuffers(int numberOfFields) {
  _numberOfFields = numberOfFields;
  if(_numberOfFields>0) {
    _id = new int[_numberOfFields];
    _area = new double[_numberOfFields];
    _x = new double[_numberOfFields];
    _y = new double[_numberOfFields];
    _e1 = new double[_numberOfFields];
    _e2 = new double[_numberOfFields];
    _pR1 = new double[_numberOfFields];
    _pR2 = new double[_numberOfFields];
    _iR = new double[_numberOfFields];
    _pG1 = new double[_numberOfFields];
    _pG2 = new double[_numberOfFields];
    _iG = new double[_numberOfFields];
    _pB1 = new double[_numberOfFields];
    _pB2 = new double[_numberOfFields];
    _iB = new double[_numberOfFields];
    _duringTransitionBefore = new int[_numberOfFields];
    _mother = new int[_numberOfFields];
    _sister = new int[_numberOfFields];
    _duringTransitionAfter = new int[_numberOfFields];
    _daughter = new int[_numberOfFields];
  }
}

Cell::DataFileIO::~DataFileIO() {
  if(_numberOfFields>0) {
    delete[] _id;
    delete[] _area;
    delete[] _x;
    delete[] _y;
    delete[] _e1;    
    delete[] _e2;    
    delete[] _pR1;    
    delete[] _pR2;    
    delete[] _iR;    
    delete[] _pG1;    
    delete[] _pG2;    
    delete[] _iG;    
    delete[] _pB1;    
    delete[] _pB2;    
    delete[] _iB;    
    delete[] _duringTransitionBefore;
    delete[] _mother;
    delete[] _sister;
    delete[] _duringTransitionAfter;
    delete[] _daughter;
  }  
}
   
bool Cell::DataFileIO::load(TissueState &s, const NcFile &dataFile) {
  // Create the cell dimension
  NcDim *cellsDim;
  if(!(cellsDim=dataFile.get_dim("cell index"))) {
    return true;
  }
  allocateBuffers(cellsDim->size());
  
  // Create the variables
  if(!(_idVar = dataFile.get_var("Cell::id"))) return false;
  if(!(_areaVar = dataFile.get_var("Cell::area"))) return false;
  if(!(_xVar = dataFile.get_var("Cell::x"))) return false;
  if(!(_yVar = dataFile.get_var("Cell::y"))) return false;
  if(!(_e1Var = dataFile.get_var("Cell::e1"))) return false;
  if(!(_e2Var = dataFile.get_var("Cell::e2"))) return false;
  if(!(_pR1Var = dataFile.get_var("Cell::pR1"))) return false;
  if(!(_pR2Var = dataFile.get_var("Cell::pR2"))) return false;
  if(!(_iRVar = dataFile.get_var("Cell::iR"))) return false;
  if(!(_pG1Var = dataFile.get_var("Cell::pG1"))) return false;
  if(!(_pG2Var = dataFile.get_var("Cell::pG2"))) return false;
  if(!(_iGVar = dataFile.get_var("Cell::iG"))) return false;
  if(!(_pB1Var = dataFile.get_var("Cell::pB1"))) return false;
  if(!(_pB2Var = dataFile.get_var("Cell::pB2"))) return false;
  if(!(_iBVar = dataFile.get_var("Cell::iB"))) return false;
  if(!(_duringTransitionBeforeVar = dataFile.get_var("Cell::duringTransitionBefore"))) return false;
  if(!(_motherVar = dataFile.get_var("Cell::mother"))) return false;
  if(!(_sisterVar = dataFile.get_var("Cell::sister"))) return false;
  if(!(_duringTransitionAfterVar = dataFile.get_var("Cell::duringTransitionAfter"))) return false;
  if(!(_daughterVar = dataFile.get_var("Cell::daughter"))) return false;

  // Create arrays and load the data into them
  if(!_idVar->get(&_id[0], _numberOfFields)) return false;
  if(!_areaVar->get(&_area[0], _numberOfFields)) return false;
  if(!_xVar->get(&_x[0], _numberOfFields)) return false;
  if(!_yVar->get(&_y[0], _numberOfFields)) return false;
  if(!_e1Var->get(&_e1[0], _numberOfFields)) return false;    
  if(!_e2Var->get(&_e2[0], _numberOfFields)) return false;    
  if(!_pR1Var->get(&_pR1[0], _numberOfFields)) return false;    
  if(!_pR2Var->get(&_pR2[0], _numberOfFields)) return false;    
  if(!_iRVar->get(&_iR[0], _numberOfFields)) return false;    
  if(!_pG1Var->get(&_pG1[0], _numberOfFields)) return false;    
  if(!_pG2Var->get(&_pG2[0], _numberOfFields)) return false;    
  if(!_iGVar->get(&_iG[0], _numberOfFields)) return false;    
  if(!_pB1Var->get(&_pB1[0], _numberOfFields)) return false;    
  if(!_pB2Var->get(&_pB2[0], _numberOfFields)) return false;    
  if(!_iBVar->get(&_iB[0], _numberOfFields)) return false;
  if(!_duringTransitionBeforeVar->get(&_duringTransitionBefore[0], _numberOfFields)) return false;
  if(!_motherVar->get(&_mother[0], _numberOfFields)) return false;
  if(!_sisterVar->get(&_sister[0], _numberOfFields)) return false;
  if(!_duringTransitionAfterVar->get(&_duringTransitionAfter[0], _numberOfFields)) return false;
  if(!_daughterVar->get(&_daughter[0], _numberOfFields)) return false;

  // Create the cells
  for(int i=0; i<_numberOfFields; ++i) {
    Cell *c = s.newCell(_id[i]);
    if(!c) {
      std::cout << "Cell::DataFileIO::load: Cell with id " << std::hex << _id[i] << " already present in network!" << std::endl;
      return false;
    }
    c->area = _area[i];
    c->r.set(_x[i], _y[i]);
    (int&)c->duringTransitionBefore = _duringTransitionBefore[i];
    c->mother = _mother[i];
    c->sister = _sister[i];
    (int&)c->duringTransitionAfter = _duringTransitionAfter[i];
    c->daughter = _daughter[i];
    c->elongation.set(_e1[i], _e2[i]);
    c->polarityB.set(_pB1[i], _pB2[i]);
    c->intIntensityB = _iB[i];
    c->polarityG.set(_pG1[i], _pG2[i]);
    c->intIntensityG = _iG[i];
    c->polarityR.set(_pR1[i], _pR2[i]);
    c->intIntensityR = _iR[i];
  }
  
  return true;  
}

bool Cell::DataFileIO::save(NcFile &dataFile, const TissueState &s) {
  // This condition is to allow for empty networks to be saved.
  // Some netCDF versions have problems with dimension length 0. They think, it's the unlimited dimension.
  // So, if there are no vertices, nothing will be written to the file. Not even the dimension.
  if(s.numberOfCells()) {
    // Create the bond dimension
    NcDim *cellsDim;
    if(!(cellsDim = dataFile.add_dim("cell index", s.numberOfCells()))) return false;

    // Create the variables
    if(!(_idVar = dataFile.add_var("Cell::id", ncInt, cellsDim))) return false;
    if(!(_areaVar = dataFile.add_var("Cell::area", ncDouble, cellsDim))) return false;
    if(!(_xVar = dataFile.add_var("Cell::x", ncDouble, cellsDim))) return false;
    if(!(_yVar = dataFile.add_var("Cell::y", ncDouble, cellsDim))) return false;
    if(!(_e1Var = dataFile.add_var("Cell::e1", ncDouble, cellsDim))) return false;
    if(!(_e2Var = dataFile.add_var("Cell::e2", ncDouble, cellsDim))) return false;
    if(!(_pR1Var = dataFile.add_var("Cell::pR1", ncDouble, cellsDim))) return false;
    if(!(_pR2Var = dataFile.add_var("Cell::pR2", ncDouble, cellsDim))) return false;
    if(!(_iRVar = dataFile.add_var("Cell::iR", ncDouble, cellsDim))) return false;
    if(!(_pG1Var = dataFile.add_var("Cell::pG1", ncDouble, cellsDim))) return false;
    if(!(_pG2Var = dataFile.add_var("Cell::pG2", ncDouble, cellsDim))) return false;
    if(!(_iGVar = dataFile.add_var("Cell::iG", ncDouble, cellsDim))) return false;
    if(!(_pB1Var = dataFile.add_var("Cell::pB1", ncDouble, cellsDim))) return false;
    if(!(_pB2Var = dataFile.add_var("Cell::pB2", ncDouble, cellsDim))) return false;
    if(!(_iBVar = dataFile.add_var("Cell::iB", ncDouble, cellsDim))) return false;
    if(!(_duringTransitionBeforeVar = dataFile.add_var("Cell::duringTransitionBefore", ncInt, cellsDim))) return false;
    if(!(_motherVar = dataFile.add_var("Cell::mother", ncInt, cellsDim))) return false;
    if(!(_sisterVar = dataFile.add_var("Cell::sister", ncInt, cellsDim))) return false;
    if(!(_duringTransitionAfterVar = dataFile.add_var("Cell::duringTransitionAfter", ncInt, cellsDim))) return false;
    if(!(_daughterVar = dataFile.add_var("Cell::daughter", ncInt, cellsDim))) return false;
    
    // fill the buffers
    allocateBuffers(s.numberOfCells());
    int i=0;
    for(TissueState::CellConstIterator it=s.beginCellIterator(); it!=s.endCellIterator(); ++it) {
      Cell *c = s.cell(it);
      _id[i] = c->id;
      _area[i] = c->area;
      _x[i] = c->r.x();
      _y[i] = c->r.y();
      _duringTransitionBefore[i] = (int)c->duringTransitionBefore;
      _mother[i] = c->mother;
      _sister[i] = c->sister;
      _duringTransitionAfter[i] = (int)c->duringTransitionAfter;
      _daughter[i] = c->daughter;
      _e1[i] = c->elongation.c1();
      _e2[i] = c->elongation.c2();
      _pB1[i] = c->polarityB.c1();
      _pB2[i] = c->polarityB.c2();
      _iB[i] = c->intIntensityB;
      _pG1[i] = c->polarityG.c1();
      _pG2[i] = c->polarityG.c2();
      _iG[i] = c->intIntensityG;
      _pR1[i] = c->polarityR.c1();
      _pR2[i] = c->polarityR.c2();
      _iR[i] = c->intIntensityR;
      ++i;
    }

    // put the records
    if(!_idVar->put(_id, s.numberOfCells())) return false;
    if(!_areaVar->put(_area, s.numberOfCells())) return false;
    if(!_xVar->put(_x, s.numberOfCells())) return false;
    if(!_yVar->put(_y, s.numberOfCells())) return false;
    if(!_e1Var->put(_e1, s.numberOfCells())) return false;
    if(!_e2Var->put(_e2, s.numberOfCells())) return false;
    if(!_pR1Var->put(_pR1, s.numberOfCells())) return false;
    if(!_pR2Var->put(_pR2, s.numberOfCells())) return false;
    if(!_iRVar->put(_iR, s.numberOfCells())) return false;
    if(!_pG1Var->put(_pG1, s.numberOfCells())) return false;
    if(!_pG2Var->put(_pG2, s.numberOfCells())) return false;
    if(!_iGVar->put(_iG, s.numberOfCells())) return false;
    if(!_pB1Var->put(_pB1, s.numberOfCells())) return false;
    if(!_pB2Var->put(_pB2, s.numberOfCells())) return false;
    if(!_iBVar->put(_iB, s.numberOfCells())) return false;
    if(!_duringTransitionBeforeVar->put(_duringTransitionBefore, s.numberOfCells())) return false;
    if(!_motherVar->put(_mother, s.numberOfCells())) return false;
    if(!_sisterVar->put(_sister, s.numberOfCells())) return false;
    if(!_duringTransitionAfterVar->put(_duringTransitionAfter, s.numberOfCells())) return false;
    if(!_daughterVar->put(_daughter, s.numberOfCells())) return false;
  }
  return true;  
}

#endif
