#ifndef BPINTERFACE_H_
#define BPINTERFACE_H_
 
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <list>

#include "boost/program_options.hpp"

#include "enum.h"
#include "vector.h"
#include "Factor.h"
#include "graphmodel.h"
#include "lbp.h"
#include "gbp.h"
#include "VarSet.h"

using mex::mxObject;
using mex::Var;
using mex::VarSet;
using mex::Factor;
using mex::vector;
using mex::graphModel;
using mex::timeSystem;
using mex::gbp;


/*
 * BPinterface.h
 * Defines methods to be called by APPRIL_Interface
 *
 *  Created on: 22 Oct 2014
 *      Author: kathryn rodgers
 */

namespace po = boost::program_options;


namespace lgbp {
#define c_log10 std::log(10)


  class BpInterface {

  protected:

  public:
    //
    // Fields
    //
    MEX_ENUM(Task , MPE,PR,MAR,MMAP );
	
    //
    // Functions
    //
    //BpInterface();
	
   
    bool initialize(int argc, char** argv); // Initialize the computation. Returns true on success.

    bool estimateComplexity(int &timeComplexity, int &memComplexity);  // Does nothing right now  
    double getPRSolution();
    mex::vector<Factor> getMARSolution();
    bool runInference(int timeLimit); // runs the inference algorithm for a small amount of time. Meant to be called repeatedly

  private:
    //
    // Fields
    //
    MEX_ENUM(Phase , EXACT, LBP, GBP, ITERCOND, DONE );

    int task; // which task we're doing
    int phase; // track if stop command has been sent
    mex::vector<Factor> bel; // stores results for MAR
    mex::vector<Factor> facts; // stores the factors 
    VarSet evVar; // evidence variables
    mex::VarOrder order; // the ordering of the variables
    size_t nvar;
    mex::graphModel factGraph; // the graph of the model
    double logZ; // result of PR

    
	
    //struct algOptions{
    double MemLimit;
    double lbpTime, lbpIter, lbpObj, lbpErr;
    double gbpTime, gbpIter, gbpObj;
    double dt;
    int    doCond;
    bool   doVerbose; // toggles printing to std::out
    int    nExtra;
    int    iboundInit;
    	
    bool isExact; // true when we've gotten an exact solution
    double timeOrder; // max time to spend on variable ordering
    int    nOrders;
    double memUseRandom;	// if memory *way* out of reach, just use random ordering...
    //} ;
	
    po::variables_map vm; // algorithm options
	
    //
    // Functions
    //
  
    bool doLoopyBP();    // Does Loopy BP on factGraph
    bool doGeneralBP();    // Does General BP on factGraph
    bool doIterativeConditioning();
    bool parseCommandOptions(int argc, char** argv);    // Puts command line options in vm
    mex::VarOrder computeVariableOrder(int numTries, double timeLimit);    // Computes a variable order
    bool fitsMBE(const graphModel& gm, const mex::VarOrder& order, const VarSet* cond);
    double solveMBE(const graphModel& gm, const mex::VarOrder& order);
    bool tryExactPR(const graphModel& gm, const mex::VarOrder& order);
    bool gbpPopulateCliques(mex::gbp& _gbp, const mex::VarOrder& order, size_t& ibound, VarSet* cond);
    bool readEvidenceFile();
    mex::vector<Factor> readUaiFile();

  };

}  // namespace lgbp

#endif /* BPINTERFACE_H_ */
