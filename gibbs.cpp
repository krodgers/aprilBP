#include <assert.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdint.h>
#include <cstdarg>
#include <cstring>

#include "factorgraph.h"

// Gibbs sampling algorithm for graphModel class
// 

using std::ifstream;
using mex::vector;
using mex::Var;
using mex::VarSet;
using mex::Factor;
using mex::graphModel;
typedef mex::graphModel::flist flist;

vector<vector<uint32_t> > gibbs(graphModel& gm, size_t nIter, size_t nSamples=1, const vector<uint32_t> init=vector<uint32_t>()) {
  vector<vector<uint32_t> > samples; samples.reserve(nSamples); // get memory for samples
  vector<uint32_t> state(init);
  if (state.size()==0) {                                        // if we need to initialize, do so
    state.resize(gm.nvar());
    for  (size_t v=0; v<gm.nvar();++v) {
      if (gm.var(v).states()==0) state[v]=0;
      else state[v] = mex::randi(gm.var(v).states());
    }
  }

  for (size_t i=0,j=0; i<nSamples; ++i) {                       // get each sample
    size_t jNext = (size_t) ((1.0 + i)/nSamples*nIter);
    for (; j<jNext; ++j) {                                      // evenly spaced over nIter
      for (size_t v=0; v<gm.nvar();++v) {                       // each iteration, go over all vars
        if (gm.var(v).states()==0) continue;                    //   (non-empty)
        const flist& factors = gm.withVariable(gm.var(v));      // get the factors they are involved with
        Factor F(gm.var(v),1.0);
        for (flist::const_iterator f=factors.begin(); f!=factors.end(); ++f) {
          VarSet vs=gm.factor(*f).vars();  vs/=gm.var(v);       // and condition on their neighbors
          //vector<uint32_t> val; val.resize(vs.size());
          //for (size_t k=0;k<vs.size();++k) val[k] = state[vs[k]];
          F *= gm.factor(*f).slice(vs, sub2ind(vs,state));
        }
        state[v] = F.sample();                                  // then draw a new value
      } 
    }
    samples.push_back(state);
  }
  return samples;
}


vector<uint32_t> annealedGibbs(graphModel& gm, size_t nIter, const vector<uint32_t> init=vector<uint32_t>()) {
  vector<uint32_t> state(init), best(init);
  double bestscore, score=0.0, temp=.5;
  bestscore = -std::numeric_limits<double>::infinity();
  if (state.size()==0) {                                        // if we need to initialize, do so
    state.resize(gm.nvar());
    for  (size_t v=0; v<gm.nvar();++v) {
      if (gm.var(v).states()==0) state[v]=0;
      else state[v] = mex::randi(gm.var(v).states());
    }
  }
  for (size_t f=0;f<gm.nFactors();++f) {
    VarSet vs=gm.factor(f).vars();
    vector<uint32_t> val; val.resize(vs.size());
    for (size_t v=0;v<vs.size();v++) val[v]=state[vs[v]];
    score += std::log( gm.factor(f)[sub2ind(vs,val)] );
  }

  for (size_t i=0,j=0; i<nIter; ++i) {                        // at each step 
    for (size_t v=0; v<gm.nvar();++v) {                       // each iteration, go over all vars
      if (gm.var(v).states()==0) continue;                    //   (non-empty)
      const flist& factors = gm.withVariable(gm.var(v));      // get the factors they are involved with
      Factor F(gm.var(v),1.0);
      for (flist::const_iterator f=factors.begin(); f!=factors.end(); ++f) {
        VarSet vs=gm.factor(*f).vars();  vs/=gm.var(v);       // and condition on their neighbors
        //vector<uint32_t> val; val.resize(vs.size());          
        //for (size_t k=0;k<vs.size();++k) val[k] = state[vs[k]];
        F *= gm.factor(*f).slice(vs, sub2ind(vs,state));
      }
      score -= std::log(F[state[v]]);
      state[v] = (F^temp).sample();                           // then draw a new value
      score += std::log(F[state[v]]);
      if (score>bestscore) { bestscore=score; best=state; }
    } 
    temp=temp*1.01;
  }
  return best;
}








#ifdef MEX
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
#else
int main(int argc, char* argv[])
#endif
{

  mex::graphModel gm;
  char fdefault[] = "test.uai";
  char* fname = fdefault;

#ifndef MEX
  if (argc>1) fname = argv[1];
#else
  std::cout<<"Pass inital factor graph if desired \n";
  std::cout<<"=========================\n\n";
  if (nrhs>0) gm.mxCopy(prhs[0]);
  else
#endif
  {
    ifstream is; is.open(fname);
    if (!is.is_open()) throw std::runtime_error("Failed to open file");
    mex::vector<Factor> flist = Factor::readUai10(is);
    for (mex::vector<Factor>::iterator i=flist.begin();i!=flist.end();++i) gm.addFactor(*i);
  }

std::cout<<"Factorgraph has "<<gm.nFactors()<<" factors over "<<gm.nvar()<<" variables\n";

for (int v=0;v<gm.nvar();++v) std::cout<<gm.var(v).label()<<":"<<gm.var(v).states()<<" ";
std::cout<<"\n";

vector<vector<uint32_t> > samples; samples = gibbs(gm,100,10);

std::cout<<samples.size()<<" samples:\n";

for (int i=0;i<samples.size();++i) {
  for (int j=0; j<samples[i].size();++j) 
    std::cout<<samples[i][j]<<" ";
  std::cout<<"\n";
}

std::cout<<"\n\n";
vector<uint32_t> best; best = annealedGibbs(gm,100);
for (int j=0; j<best.size();++j) std::cout<<best[j]<<" "; std::cout<<"\n";


std::cout<<"=========================\n\n\n";

}


