/////////////////////////////////////////////////////////////////////////////////////
// Factor.h  --  class definition for matlab-compatible factor class
//
// A few functions are defined only for MEX calls (construction & load from matlab)
// Most others can be used more generally.
// 
//////////////////////////////////////////////////////////////////////////////////////


#ifdef WINDOWS
//#include "stdafx.h"
#include <iostream>


#define strcasecmp _stricmp
#endif

//#include <assert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <algorithm>
#include <numeric>
#include <float.h>
#include <vector>

#include "mxUtil.h"
#include "VarSet.h"
#include "subindex.h"
#include "vector.h"

#include "Factor.h"

namespace mex {


#ifdef __FACTOR_H_MEMORY
size_t Factor::memused = 0;
size_t Factor::mmax = 0;
#endif

/************************************************************************************
 ************************************************************************************
 IMPLEMENTATION
 ************************************************************************************
************************************************************************************/

//////////////////////////////////////////////////////////////////////////////////////////////
// MEX specific functions
//////////////////////////////////////////////////////////////////////////////////////////////
#ifdef MEX
bool Factor::mxCheckValid(const mxArray* M) {
  return mxIsStruct(M) && t_.mxCheckValid(mxGetField(M,0,"t")) && v_.mxCheckValid(mxGetField(M,0,"v"));
}

void Factor::mxSet(mxArray* M) {
  if (mxAvail()) {
    //mxDestroyArray(M_);            //   destroy old matlab object if it existed
  } // !!! what to do about v_ and t_?
  M_=M;
  if (M) {                            // if we got a matlab object, set to it
    t_.mxSetND(mxGetField(M_,0,"t"));
    v_.mxSet(mxGetField(M_,0,"v"), mxGetDimensions(t_.mxGet()));
  } else {                            // null pointer => clear 
    *this = Factor();
  }
}

mxArray* Factor::mxGet() {
  if (!mxAvail()) {
    mxArray* m; 
    int retval = mexCallMATLAB(1,&m,0,NULL,"factor");           // create new empty factor object
    if (retval) mexErrMsgTxt("Error creating new factor");      //   associated with a matlab mxArray
    Factor f; f.mxSet(m);                                        // use the copy constructor to put our data
    f = *this;                                                  //   there (in matlab-allocated memory)
    mxSwap(f);                                                  // then swap everything for the new object
    setDims();
  }
  return M_;
};

void Factor::mxRelease() {                                 // Disassociate with a given matlab object
  if (mxAvail()) {                                         //   without deallocating / destroying it
    t_.mxRelease(); v_.mxRelease(); 
  }
}
void Factor::mxDestroy() {                                  //  Disassociate with a given matlab object
  mxArray* m=M_;                                            //    and also destroy it
  mxRelease();
  if (m) mxDestroyArray(m);
}
void Factor::mxSwap(Factor& f) {
  v_.mxSwap(f.v_);
  t_.mxSwap(f.t_);
  std::swap(M_,f.M_);
}
#endif
//////////////////////////////////////////////////////////////////////////////////////////////


vector<Factor> Factor::readUai10(std::istream& is) {
  size_t nvar, ncliques, csize, v, nval;
  char* st; st = new char[20];
  is >> st;
 
  if ( strcasecmp(st,"MARKOV") ) 
    throw std::runtime_error("Only UAI-2010 Markov-format files are supported currently");

  is >> nvar;
  vector<size_t> dims(nvar);
  for (size_t i=0;i<nvar;i++) { is>>dims[i]; if (dims[i] == 0) dims[i]=1; }
  
  is >> ncliques;
  std::vector<mex::vector<Var> > cliques(ncliques);
  std::vector<VarSet > sets(ncliques);
  for (size_t i=0;i<ncliques;i++) {
    is >> csize;
    cliques[i].reserve(csize);
    for (size_t j=0;j<csize;j++) { 
      is>>v; Var V(v,dims[v]);
      cliques[i].push_back(V); 
      sets[i] |= V;
    }   
  }
  
  //vector<vector<double> > tables(ncliques);
  vector<Factor> tables(ncliques);
  for (size_t i=0;i<ncliques;i++) {
    is >> nval;
    assert(nval == sets[i].nrStates());
    tables[i] = Factor(sets[i],0.0);        // preallocate memory and convert from given order, bigEndian
		if (nval==1) { is>>tables[i][0]; continue; }  // special case for constant factor
    permuteIndex pi( cliques[i], true); pi=pi.inverse();   // to our order
    for (size_t j=0;j<nval;j++) is>>tables[i][pi.convert(j)];
  }

  delete[] st;
  return tables;
}


void Factor::writeUai10(std::ostream& os, const vector<Factor>& fs) {
  size_t nvar=0;

  for (size_t f=0;f<fs.size();++f) if (fs[f].nvar()>0)          // find maximum variable label
    nvar=std::max(nvar,(size_t)(fs[f].vars().rbegin()->label()+1));
  
  vector<uint32_t> dims(nvar,1);                                // collect all variable dimensions
  for (size_t f=0;f<fs.size();++f) 
    for (size_t v=0;v<fs[f].nvar();++v) {
      Var V=fs[f].vars()[v]; dims[V.label()] = V.states(); 
    }

  os << "MARKOV\n";                                             // Markov Random Field:
  os << nvar << "\n";                                           // # variables
  for (size_t v=0;v<nvar;++v) os << dims[v] << " ";  os<<"\n";  // dimensions
  os << fs.size() << "\n";                                      // # cliques
  for (size_t f=0;f<fs.size();++f) {
    os << fs[f].nvar() << "   ";                                // each clique: # var, var ids
    for (size_t v=0;v<fs[f].nvar();++v) os << fs[f].vars()[v].label() << " ";
    os << "\n";
  }
  for (size_t f=0;f<fs.size();++f) {                            // each function
    os << fs[f].nrStates() << "   ";                            // # states, then values
		if (fs[f].nrStates()==1) { os <<std::scientific << fs[f][0] <<"\n"; continue; } // special case: no vars
    vector<Var> clique( fs[f].vars().begin(), fs[f].vars().end() );
    permuteIndex pi( clique, true); pi=pi.inverse();
    for (size_t j=0;j<fs[f].nrStates();++j) os << std::scientific << fs[f][pi.convert(j)] << " ";
    os << "\n";
  }
}



} // end namespace mex

