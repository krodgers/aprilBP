/////////////////////////////////////////////////////////////////////////////////////
// Factor.h  --  class definition for matlab-compatible factor class
//
// A few functions are defined only for MEX calls (construction & load from matlab)
// Most others can be used more generally.
// 
//////////////////////////////////////////////////////////////////////////////////////
//
// Written by Alex Ihler
// Copyright (C) 2010 Alexander Ihler; distributable under GPL -- see README.txt
//
//////////////////////////////////////////////////////////////////////////////////////
#ifndef __FACTOR_H
#define __FACTOR_H

// MEMORY define: keep track of global Factor memory usage (helpful in debugging other code)
//#define __FACTOR_H_MEMORY

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
#include "enum.h"
#include "VarSet.h"
#include "subindex.h"
#include "vector.h"

namespace mex {

class Factor : public virtual mxObject {
 public:
  typedef double value;
  typedef VarSet::vindex vindex;   // variable identifiers (1...N)
  typedef VarSet::vsize  vsize;    // variable values (0...K-1)

  // Constructors //////////////////////////////////////////////////////////////////////////

  Factor(Factor const& f) : v_(f.v_), t_(f.t_) {
#ifdef __FACTOR_H_MEMORY
    memused += t_.capacity()*sizeof(double);
    mmax = std::max(mmax,memused);
#endif
  }         // copy ctor
  explicit Factor(const value s=1.0) : v_(),   t_(1,s)  { 
#ifdef __FACTOR_H_MEMORY
    memused += t_.capacity()*sizeof(double);
    mmax = std::max(mmax,memused);
#endif
  }        // scalar constructor
  explicit Factor(VarSet const& vs, value s=1.0) : v_(vs),t_() {   // constant factor over given vars
    t_.resize(vs.nrStates()); fill(s);
#ifdef __FACTOR_H_MEMORY
    memused += t_.capacity()*sizeof(double);
    mmax = std::max(mmax,memused);
#endif
  }
  Factor(VarSet const& vs, value* T) : v_(vs),t_() { 
    t_.resize(v_.nrStates()); std::copy(T,T+t_.size(),t_.begin()); 
#ifdef __FACTOR_H_MEMORY
    memused += t_.capacity()*sizeof(double);
    mmax = std::max(mmax,memused);
#endif
  }
  // Factor( from Vars and a vector of values)  !!
  // Factor( from Vars and a value* )           !!
  // Factor( from permuted vector of Var and vector of values) !!

  ~Factor() { 
#ifdef __FACTOR_H_MEMORY
    memused -= t_.capacity()*sizeof(double); 
#endif
   }                                    // destructor

  // Assignments & copy constructors //////////////////////////////////////////////////////
  
  Factor& operator=(Factor const& rhs) {                  // assignment (deep copy)
    if (this!=&rhs) { 
#ifdef __FACTOR_H_MEMORY
      memused -= t_.capacity()*sizeof(double);
#endif
#ifndef MEX
      { vector<double> tmp; t_.swap(tmp); }                // force vector to release memory
#endif
      v_ = rhs.v_; t_ = rhs.t_; setDims();                 // then reassign
#ifdef __FACTOR_H_MEMORY
      memused += t_.capacity()*sizeof(double);
      mmax = std::max(mmax,memused);
#endif
    }
    return *this;
  }

#ifdef __FACTOR_H_MEMORY
  static size_t memused;
  static size_t mmax;
  static size_t mem() { return memused; };
#endif

  void swap(Factor& F) {                                  // object contents exchange
    if (&F!=this) { v_.swap(F.v_); t_.swap(F.t_); setDims(); F.setDims(); }
  }

  virtual Factor* clone() const { Factor *F=new Factor(*this); return F; }  // clone from pointer
  
  // MEX Class Wrapper Functions //////////////////////////////////////////////////////////
#ifdef MEX 
   bool        mxCheckValid(const mxArray*);   // check if matlab object is compatible with this object
   void        mxSet(mxArray*);     // associate with A by reference to data
   mxArray*    mxGet();             // get a pointer to the matlab object wrapper (creating if required)
   void        mxRelease();         // disassociate with a matlab object wrapper, if we have one
   void        mxDestroy();         // disassociate and delete matlab object
   void        mxSwap(Factor&);     // disassociate and delete matlab object
#endif
  void setDims() {                  // make variable & table dimensions consistent (matlab only)
#ifdef MEX
    if (t_.mxAvail())
      if (v_.size()>1) mxSetDimensions(t_.mxGet(),v_.dims(),v_.size());
      else { mxSetM(t_.mxGet(),v_.size() ? v_.dims()[0] : 1); mxSetN(t_.mxGet(),1); }
#endif
  }

  /////////////////////////////////////////////////////////////////////////////////////////
  // Accessor Functions ///////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////
  const size_t     nvar()      const { return v_.nvar();  };// # of variables
  const VarSet&    vars()      const { return v_;   };      // variable IDs
  // Non-const version removed; unsafe? !!
  const VarSet&    variables() const { return v_;   };      // name difference (!!)
  const vsize*     dims()      const { return v_.dims();};  // variable dimensions
  const value*     table()     const { return &t_[0];   };  // table of factor values
  size_t           nrStates()  const { return t_.size(); }  // get table size
  size_t           numel()     const { return t_.size(); }; // name difference (!!)

  // Boolean checks on factor properties //////////////////////////////////////////////////
  bool isempty()  const { return t_.empty(); };             // empty factor
  //bool isnan()    const { return std::find_if(t_.begin(),t_.end(); isnan) != t_.end(); }; // !!!
  bool isnan()    const { bool b=false; for (size_t i=0;i<nrStates();i++) b|=isnan(t_[i]);    return b; };
  bool isfinite() const { bool b=true;  for (size_t i=0;i<nrStates();i++) b&=isfinite(t_[i]); return b; };
  bool isscalar() const { return numel()==1; };

  void cleanInf() { for (size_t i=0;i<nrStates();i++) if (!isfinite(t_[i])) t_[i]=0.0; };

  // Direct value accessor ////////////////////////////////////////////////////////////////
  value  operator[] (vsize v) const { return t_[v]; };  // access table elements
  value& operator[] (vsize v)       { return t_[v]; };  // rewrite table elements
  value  get(vsize i)         const { return t_.at(i); }  // "safe" versions
  void   set(vsize i, value v)      { t_.at(i)=v; }

  // Filling
  Factor& fill(value v) { std::fill(t_.begin(),t_.end(),v); return *this; }  // fill with constant
  Factor& randomize()   { std::generate(t_.begin(),t_.end(),mex::randu); return *this; }
  Factor& fillUniform() { return fill(1.0/nrStates()); }                  // fill with 1/N

  // Unary transformations (in-place); see outside class def'n for const versions /////////
  inline Factor& abs(void)    { std::transform(t_.begin(),t_.end(),t_.begin(),unOpAbs()); return *this; };
  inline Factor& exp(void)    { std::transform(t_.begin(),t_.end(),t_.begin(),unOpExp()); return *this; };
  inline Factor& log0(void)   { std::transform(t_.begin(),t_.end(),t_.begin(),unOpLog0()); return *this; };
  inline Factor& log(void)    { std::transform(t_.begin(),t_.end(),t_.begin(),unOpLog()); return *this; };
  inline Factor& log2(void)   { std::transform(t_.begin(),t_.end(),t_.begin(),unOpLogL(std::log(2))); return *this; };
  inline Factor& log10(void)  { std::transform(t_.begin(),t_.end(),t_.begin(),unOpLog10()); return *this; };
  // Normalize? (!!)

  // Functors defined for unary operations (transformations) to the factor table
  struct unOpAbs   { value operator()(value a)  { return std::abs(a);  } };
  struct unOpExp   { value operator()(value a)  { return std::exp(a);  } };
  struct unOpLog0  { value operator()(value a)  { return (a!=0.0) ? std::log(a) : 0.0;  } };
  struct unOpLog   { value operator()(value a)  { return std::log(a);  } };
  struct unOpLogL  { value l; unOpLogL(value L):l(L) {}; value operator()(value a)  { return std::log(a)/l;  } };
  struct unOpLog10 { value operator()(value a)  { return std::log10(a);} };
  
  struct myMin { value operator()(value a,value b)  { return std::min(a,b);} };
  struct myMax { value operator()(value a,value b)  { return std::max(a,b);} };

////////////////////////////////////////////////////////////////////////////////
// Basic factor operations (+,-,*,/)
////////////////////////////////////////////////////////////////////////////////
  Factor  operator+ (const Factor& B) const  { return binaryOp(  B, binOpPlus()  ); };
  Factor& operator+=(const Factor& B)        { return binaryOpIP(B, binOpPlus()  ); };
  Factor  operator- (const Factor& B) const  { return binaryOp(  B, binOpMinus() ); };
  Factor& operator-=(const Factor& B)        { return binaryOpIP(B, binOpMinus() ); };
  Factor  operator* (const Factor& B) const  { return binaryOp(  B, binOpTimes() ); };
  Factor& operator*=(const Factor& B)        { return binaryOpIP(B, binOpTimes() ); };
  Factor  operator/ (const Factor& B) const  { return binaryOp(  B, binOpDivide()); };
  Factor& operator/=(const Factor& B)        { return binaryOpIP(B, binOpDivide()); };

  Factor  operator+ (const value B) const    { return binaryOp(  B, binOpPlus()  ); };
  Factor& operator+=(const value B)          { return binaryOpIP(B, binOpPlus()  ); };
  Factor  operator- (const value B) const    { return binaryOp(  B, binOpMinus() ); };
  Factor& operator-=(const value B)          { return binaryOpIP(B, binOpMinus() ); };
  Factor  operator* (const value B) const    { return binaryOp(  B, binOpTimes() ); };
  Factor& operator*=(const value B)          { return binaryOpIP(B, binOpTimes() ); };
  Factor  operator/ (const value B) const    { return binaryOp(  B, binOpDivide()); };
  Factor& operator/=(const value B)          { return binaryOpIP(B, binOpDivide()); };
  Factor  operator^ (const value B) const    { return binaryOp(  B, binOpPower() ); };
  Factor& operator^=(const value B)          { return binaryOpIP(B, binOpPower() ); };

  // Above operators use the following internal definitions:
  
  // Binary operations (eg A + B); returns new object
  template<typename Function> Factor binaryOp( const Factor& B, Function Op) const {
    VarSet v = v_ + B.v_;              // expand scope to union
    Factor F(v);                         //  and create target factor
    subindex s1(v,v_), s2(v,B.v_);           // index over A and B & do the op
    //for (vector<value>::iterator i=F.t_.begin(); i!=F.t_.end(); ++i,++s1,++s2) *i=Op(t_[s1], B[s2]);
    for (size_t i=0; i<F.nrStates(); ++i,++s1,++s2) F[i]=Op(t_[s1], B[s2]);
    return F;                     // return the new copy
  }
  template<typename Function> Factor binaryOp( const value B, Function Op) const {
    Factor F=*this; F.binaryOpIP(B,Op); return F; // for scalar args, define with an in-place operator
  }
  // Binary in-place operations (eg A += B); returns reference to modified A
  template<typename Function> Factor& binaryOpIP( const Factor& B, Function Op) {
    if (!(v_ >> B.v_)) { Factor F=binaryOp(B,Op); swap(F); }  // if A's scope is too small, call non-in-place version
    //VarSet v = v_ + B.v_;                  // expand scope to union
    //if (v != v_) *this = binaryOp(B,Op);           // if A's scope is too small, call binary op
    else {                         
      subindex s2(v_,B.v_);                       // otherwise create index over B
      //for (vector<value>::iterator i=t_.begin(); i!=t_.end(); ++i,++s2) *i=Op.IP(*i, B[s2]);
      for (size_t i=0; i<nrStates(); ++i,++s2) Op.IP(t_[i] , B[s2]);  // and do the operations
    }                                               
    return *this; 
  }
  template<typename Function> Factor& binaryOpIP( const value B, Function Op) {
    for (size_t i=0;i<nrStates();i++) Op.IP(t_[i] , B); return *this;  // simplifies for scalar args
  }

  // Functors defined for binary operations on the factor table : Op(a,b) and Op.IP(a,b) (in-place version)
  struct binOpPlus   { 
    value  operator()(value  a, const value b) { return a+b; }; 
    value&         IP(value& a, const value b) { return a+=b;}; 
  };
  struct binOpMinus  { 
    //value  operator()(value  a, const value b) { return a-b; }; 
    //value&         IP(value& a, const value b) { return a-=b;}; 
    value  operator()(value  a, const value b) { return (b!=-infty()) ? a-b : b; }; 
    value&         IP(value& a, const value b) { return (b!=-infty()) ? a-=b: a=b; }; 
    //value  operator()(value  a, const value b) { return (isnan(a-=b)) ? a=-infty() : a; };
    //value&         IP(value& a, const value b) { return (isnan(a-=b)) ? a=-infty() : a; };
    //value  operator()(value  a, const value b) { value c=a-b; if (c!=c) c=a; return c; }; 
    //value&         IP(value& a, const value b) { value c=a-b; if (c==c) a=c; return a; };
  };
  struct binOpTimes  { 
    value  operator()(value  a, const value b) { return a*b; }; 
    value&         IP(value& a, const value b) { return a*=b;}; 
  };
  struct binOpDivide { 
    value  operator()(value  a, const value b) { return (b) ? a/b  : 0;  }; 
    value&         IP(value& a, const value b) { return (b) ? a/=b : a=0;}; 
  };
  struct binOpPower { 
    value  operator()(value  a, const value b) { return std::pow(a,b);  }; 
    value&         IP(value& a, const value b) { return a=std::pow(a,b);}; 
  };


////////////////////////////////////////////////////////////////////////////////
// Comparative thresholding operators (useful? !!!)
////////////////////////////////////////////////////////////////////////////////
// operator> (const Factor& F)
// operator>=(const Factor& F)
// operator< (const Factor& F)
// operator<=(const Factor& F)


////////////////////////////////////////////////////////////////////////////////
// Partition function, entropy, and normalization
////////////////////////////////////////////////////////////////////////////////
  Factor  normalized()  const { Factor F=*this; F.normalize(); return F;  }
  Factor& normalize() { double Z=sum(); if (Z!=0) *this/=Z; return *this; }

  double logpartition() const { return std::log(sum()); };
  //double logpartition() const { return std::log(sum())/std::log(2.0); };

  value entropy(void) const { 
    value H=0, Z=0;
    for (size_t i=0;i<nrStates();i++) { 
      Z += t_[i]; 
      double L=std::log(t_[i]);
      if (!std::isinf(L)) H -= t_[i]*L; //std::log(t_[i]); 
      //!!!! if (!isinf(L)) H -= t_[i]*L; //std::log(t_[i]); 
    }
    H/=Z; H+=std::log(Z);
    return H;
  }


////////////////////////////////////////////////////////////////////////////////
// Elimination operators (sum, max, min, ...)
////////////////////////////////////////////////////////////////////////////////

  Factor sum(VarSet const& sumOut) const { VarSet t=v_ - sumOut; return marginal(t); }
  value  sum()                     const { return std::accumulate(t_.begin(),t_.end(),0.0,std::plus<value>()); }

  Factor sumPower(VarSet const& sumOut,value pow)  const {  
    if (pow==1.0)          return sum(sumOut); 
    else if(pow==-infty()) return min(sumOut);
    else if(pow== infty()) return max(sumOut); 
    else {
      Factor F=*this; F.log(); F*=pow; F=F.logsumexp(sumOut); F/=pow; F.exp();
      return F;
    }
  }

  Factor logsumexp(const VarSet& sumOut) const {  
    VarSet target = v_ - sumOut;
    Factor mx = maxmarginal(target);
    Factor B(target,0.0);
    subindex s2(v_,B.v_);		// !!! TODO: change to more memory efficient superindex version? see gbp.h
    for (size_t i=0; i<nrStates(); ++i,++s2) if (mx[s2]!=-infty()) B[s2]+=std::exp(t_[i] - mx[s2]);
    for (size_t i=0; i<B.nrStates();++i) mx[i] += std::log(B[i]);
    return mx;
  }
  double logsumexp() const { 
    double r=0, mx=max(); 
    if (mx == -infty()) return mx;
    for (size_t i=0;i<nrStates();++i) r+=std::exp(t_[i]-mx);
    return std::log(r)+mx;
  }

/*  // Alternative version (?)
  Factor logsumexp(const VarSet& sumOut) const {  
  VarSet target = v_ - sumOut;
  Factor F(target,0.0);
  subindex s(v_,target);
  for (size_t i=0;i<nrStates(); ++i,++s) { 
    value mx=F[s],mn=t_[i];  if (mx<mn) { mx=t_[i]; mn=F[s]; }
    F[s]=mx+std::log(1+std::exp(mn-mx)); 
  }
  return F;
  }
*/

  Factor max(VarSet const& sumOut) const { VarSet t=v_ - sumOut; return maxmarginal(t); }
  value  max()                     const { return std::accumulate(t_.begin(),t_.end(),-infty(),myMax()); }

  Factor min(VarSet const& sumOut) const { VarSet t= v_ - sumOut; return minmarginal(t); }
  value  min()                     const { return std::accumulate(t_.begin(),t_.end(),-infty(),myMin()); }
  
  size_t argmax() const { return std::distance(t_.begin(), std::max_element(t_.begin(),t_.end())); }
  size_t argmin() const { return std::distance(t_.begin(), std::min_element(t_.begin(),t_.end())); }

  size_t argmax(const VarSet& vCond, vsize vState) const { 
    if (vCond.size()==0) return argmax();
    VarSet vKeep = vars()-vCond;
		superindex sup(vars(),vKeep,vState); size_t N=vKeep.nrStates();
    size_t mxi=0; double mx=-infty();
    for (size_t i=0;i<N;++i,++sup) if (t_[sup] > mx) { mx=t_[sup]; mxi=sup; }
		return mxi;

    //superindex sup(vars(),vars()-vCond,vState);
    //size_t mxi=sup; double mx=t_[sup];
    //for (size_t i=0;sup!=sup.end();++i,++sup) if (t_[sup] > mx) { mx=t_[sup]; mxi=sup; } 
    ////ALT: for (size_t i=0;i<vKeep.nrStates() && sup!=sup.end();++i,++sup) { std::cout<<sup<<" "; if (t_[sup] > mx) { mx=t_[sup]; mxi=sup; } }

    //OLD version:
    //subindex src(vars(),vCond);
    //size_t mxi=0; double mx=-infty();
    //for (size_t i=0;i<nrStates();++i,++src) if (src==vState) if (t_[i] > mx) { mx=t_[i]; mxi=i; }
    //return mxi;
  }

  Factor condition(const VarSet& vRem, vsize vState) const { 
    VarSet vKeep = vars() - vRem;
    Factor F(vKeep,0.0);
    superindex sup(vars(),vKeep,vState); size_t N=vKeep.nrStates();
    for (size_t i=0;i<N;++i,++sup) F[i]=t_[sup];
    //for (size_t i=0;sup!=sup.end();++i,++sup) F[i]=t_[sup];
		// old version:
    //subindex src(vars(),vRem), dst(vars(),vKeep); 
    //for (size_t i=0;i<nrStates();++i,++src,++dst) if (src==vState) F[dst]=t_[i];  // !!! terrible; needs supindex
    return F;
  }
  Factor slice(const VarSet& vRem, vsize vState) const { return condition(vRem,vState); } // !! name change

  Factor embed(const VarSet& v) const { if (vars()==v) return *this; else return (*this)+Factor(v/vars(),0.0);}

  size_t sample() const { 
    double x=0.0, z = sum();  
    if (z==0.0) return mex::randi(nrStates()); 
//{ size_t r = mex::randi(nrStates()); std::cout<<r; return r; } // return (nrStates()*mex::randu()); 
		double y = mex::randu() * z;
    for (size_t i=0;i<nrStates();++i) 
      if ((x+=t_[i]) > y) return i;
    return nrStates()-1;
  }
  vsize draw() const { return sample(); }  // !! name difference

  Factor marginal(VarSet const& target) const {  
    Factor F(target&vars(),0.0);
    subindex s(v_,F.vars());
    for (size_t i=0;i<nrStates(); ++i,++s) F[s]+=t_[i];
    return F;
  }
  void marginalInto(VarSet const& target, Factor& F) const {  
    assert(F.vars()==(target&vars()) && "marginalInto: target factor has incorrect variables");
    F.fill(0.0);
    subindex s(v_,F.vars());
    for (size_t i=0;i<nrStates(); ++i,++s) F[s]+=t_[i];
    //return F;
  }

  Factor maxmarginal(VarSet const& target) const {
    Factor F(target&vars(),-infty());
    subindex s(v_,F.vars());
    for (size_t i=0;i<nrStates(); ++i,++s) F[s]=(F[s] < t_[i]) ? t_[i] : F[s];
    return F;
  }

  Factor minmarginal(VarSet const& target) const {
    Factor F(target&vars(),infty());
    subindex s(v_,F.vars());
    for (size_t i=0;i<nrStates(); ++i,++s) F[s]=(F[s] > t_[i]) ? t_[i] : F[s];
    return F;
  }


////////////////////////////////////////////////////////////////////////////////
// Misc other functions
////////////////////////////////////////////////////////////////////////////////

  MEX_ENUM( Distance , L1,L2,LInf,KL,HPM,MAS,OptGap );
  /*
  struct Distance {
    enum Type { L1, L2, LInf, KL, HPM, MAS };
    Type t_;
    Distance(Type t=L2) : t_(t) {}
    operator Type () const { return t_; }
    operator char const* () const { return names()[t_]; }
  private:
    template<typename T> operator T() const;
    static char const* const* names() {
      static char const* const str[] = {"L1","L2","LInf","KL","HPM","MAS",0};
      return str;
    }
  };
  */

  double distance(Factor const& F2, Distance type=Distance::L2) const {
    assert( vars() == F2.vars() && "distance: factors scopes do not match");
    Factor F(*this),Ftmp;               // make a copy for manipulation
    double dist=-1.0;                   // local variables
    value Z;
    switch (type) {
      case Distance::L2:                // L2, sum of squared errors
        F-=F2; F*=F; dist=F.sum();
        break;
      case Distance::L1:                // L1, sum of absolute errors  (=TV !!)
        F-=F2; dist=F.abs().sum();
        break;
      case Distance::LInf:              // L-infinity, max absolute error
        F-=F2; dist=F.abs().max();
        break;
      case Distance::KL:                // KL-divergence (relative entropy)
        Z=sum(); F/=F2; F*=F2.sum()/Z; F.log(); F*=*this; dist=F.sum()/Z;
        break;
      case Distance::HPM:               // Hilbert's projective metric
        F/=F2; F.log(); dist=F.max()-F.min(); //   aka "dynamic range"
        break;
      case Distance::MAS:               // "MAS" error value (not a metric)
        F.log(); Ftmp=F2; F/=Ftmp.log();
        dist = std::max( F.max(), 1.0/F.min() )-1.0;
        break;
      case Distance::OptGap:            // "Primal/Dual Gap"-like
        double mx1,mx2,gap1,gap2;
        mx1=mx2=gap1=gap2=0.0;
        for (size_t i=0;i<nrStates();++i) {
          if (mx1<F[i]) { gap2=F2[i]; mx1=F[i];} else if (mx1==F[i]) gap2 = std::min(gap2,F2[i]);
          if (mx2<F2[i]){ gap1=F[i]; mx2=F2[i];} else if (mx2==F2[i]) gap1= std::min(gap1,F[i]); 
        }
        return (mx1-gap1)+(mx2-gap2);
        break;
      //case Distance::Hellinger:   (!!)
      //  F^=0.5; F-=F2^0.5; F*=F; dist=(0.5*F.sum())^0.5;  // straightforward computation
      //  F*=F2; F^=0.5; dist=(1-F.sum())^0.5;              // alternate computation if F,F2 normalized
      //  break;
      default: throw std::runtime_error("Invalid distance type");
    }
    return dist;
  }

  double norm(Distance type=Distance::L2) const {
  Factor F(*this);                // make a copy for manipulation
  double dist=-1.0;               // 
  switch (type) {
    case Distance::L2:                // L2, sum of squared errors
      F*=F; dist=F.sum();
      break;
    case Distance::L1:                // L1, sum of absolute errors
      dist=F.abs().sum();
      break;
    case Distance::LInf:                // L-infinity, max absolute error
      dist=F.abs().max();
      break;
    case Distance::KL:                // KL-divergence (relative entropy => entropy?)
      return entropy();
      break;
    case Distance::HPM:               // Hilbert's projective metric
      F.log(); dist=F.max()-F.min();      //   aka "dynamic range"
      break;
    default:
      throw std::runtime_error("Invalid norm type");
  }
  return dist;
  }

  MEX_ENUM( Decomp , L2,L2_HPM,L2_MAS );

  std::vector<Factor> decompSum(std::vector<VarSet> vlist, Factor::Decomp method) const {
    int nF=vlist.size();
    double mx,mn;
    std::vector<Factor> Flist(nF);

    Factor tmp,F=*this;
    switch (method) {
      case Decomp::L2: //L2
        double Cn,Cd;
        Cd=F.numel(); Cn=F.sum(); // /Cd*(1-1.0/nF);
        for (int j=0;j<nF;j++) {
          Flist[j] = F.marginal( vlist[j] );
          double D = Cd/Flist[j].numel();
          Flist[j]/= D;
          Flist[j]-= Cn/Cd*(1.0-1.0/(nF-j));
          F  -= Flist[j];
          Cn -= Flist[j].sum()*D;
        }
        break;
      case Decomp::L2_HPM: //L2+HPM
        Flist = decompSum(vlist,Decomp::L2);
        for (int j=0;j<nF;j++) F-=Flist[j];
        mx=F.max(); mn=F.min();
        for (int j=0;j<nF;j++) Flist[j]+=(mx+mn)/2/nF;
        break;
      case Decomp::L2_MAS: //L2+MAS
        Flist = decompSum(vlist,Decomp::L2);
        F=Flist[0]; for (int j=1;j<nF;j++) F+=Flist[j];
        F /= *this; F.log();
        mx=F.max(); mn=F.min();
        for (int j=0;j<nF;j++) Flist[j]*=std::exp(-(mx+mn)/2/nF);
        break;
      }
    return Flist;
  }

  std::vector<Factor> decompProd(std::vector<VarSet> vlist, Factor::Decomp method) const {
    Factor F=*this; F.log();
    std::vector<Factor> Flist = F.decompSum(vlist,method);
    for (size_t j=0;j<vlist.size();j++) Flist[j].exp();
    return Flist;
  }
 
  friend std::ostream& operator<<( std::ostream& out, const mex::Factor& F) { 
    out << "Factor over " << F.variables() << ":";
    for (size_t j=0;j<F.t_.size();j++) out<<" "<<F.t_[j];
    return out;
  };


  static vector<Factor> readUai10(std::istream& is);
  static void           writeUai10(std::ostream& os, const vector<Factor>&);
  static vector<Factor> readErgo(std::istream& is);
  static void           writeErgo(std::ostream& os, const vector<Factor>&);
  static vector<Factor> readWCSP(std::istream& is);
  static void           writeWCSP(std::ostream& os, const vector<Factor>&);


  /////////////////////////////
  // Private object functions
  /////////////////////////////
 protected:
  VarSet v_;          // variable list vector
  vector<value> t_;    // table of values

  //vsize calcNumEl() const { size_t n=std::accumulate(dims(),dims()+nvar(),(size_t)1,std::multiplies<size_t>()); return (n>1)?n:1; }
  vsize calcNumEl() const {vsize n=1; vsize const* d=dims(); for (size_t i=0;i<nvar();i++) n*=d[i]; return (n>1)?n:1; }
#ifdef MEX 
  static inline bool isfinite(value v) { return mxIsFinite(v); }
  static inline bool isnan(value v)    { return mxIsNaN(v); }
  static inline value infty()          { return mxGetInf(); }
#else
  //static inline bool isfinite(value v) { return std::abs(v)!=std::numeric_limits<value>::infinity(); }
  static inline bool isfinite(double v) { return (v <= DBL_MAX && v >= -DBL_MAX); }
  static inline bool isnan(value v)    { return (v!=v); }
  static inline value infty()          { return std::numeric_limits<value>::infinity(); }
#endif

 public:
  static inline Factor delta(const VarSet& vs, size_t idx) { Factor F(vs,0.0); F[idx]=1.0; return F; }
};

////////////////////////////////////////////////////////////////////////////////
// "Static" functions that operate on Factor class variables
////////////////////////////////////////////////////////////////////////////////

inline Factor  operator+ (const Factor::value B, const Factor& A) { return A+B; };
inline Factor  operator* (const Factor::value B, const Factor& A) { return A*B; };
inline Factor  operator- (const Factor::value B, const Factor& A) { Factor BF(A.vars(),B); return BF-=A; };
inline Factor  operator/ (const Factor::value B, const Factor& A) { Factor BF(A.vars(),B); return BF/=A; };

inline Factor abs(const Factor& A)   { Factor F=A; F.abs(); return F; };
inline Factor exp(const Factor& A)   { Factor F=A; F.exp(); return F; };
inline Factor log0(const Factor& A)  { Factor F=A; F.log0(); return F; };
inline Factor log(const Factor& A)   { Factor F=A; F.log(); return F; };
inline Factor log2(const Factor& A)  { Factor F=A; F.log(); F/=std::log(2.0); return F; };
inline Factor log10(const Factor& A) { Factor F=A; F.log10(); return F; };

template <class InputIterator>
inline Factor mean(InputIterator first, InputIterator last) { size_t N=0; Factor F(0.0); for (;first!=last;++first,++N) F+=*first; if (N) F/=N; return F; }
template <class InputIterator>
inline Factor geomean(InputIterator first, InputIterator last) { size_t N=0; Factor F(1.0); for (;first!=last;++first,++N) F*=*first; if (N) F^=1.0/N; return F; }



} // end namespace mex

#endif
