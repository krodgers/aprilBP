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

    typedef struct {
      double MemLimit;
      double lbpTime, lbpIter, lbpObj, lbpErr;
      double gbpTime, gbpIter, gbpObj;
      int    doCond;
      bool   doVerbose; // toggles printing to std::out
      int    iboundInit; // initial ibound liimit
      double timeOrder; // max time to spend on variable ordering
      int    nOrders; // max number of variable orderings to try
      std::string problemFile; // the name of the file to read the problem from
      int task; // which task we're doing
      std::string orderFile; // name of the file to read a variable ordering from
      std::string evidenceFile; // name of the file to read evidence from
      std::string query; // name of the file to read MMAP query from
      int  nExtra;
    } algOptions; 

  class BpInterface {

  protected:

  public:
    //
    // Fields
    //
    MEX_ENUM(Task, PR, MAR );
	
    //
    // Functions
    //
    //BpInterface();
	
   
    bool initialize(int argc, char** argv); // Initialize the computation. Returns true on success.
    bool initialize(algOptions opts, bool useDefault, double totalTime); // Initialize with given paramters; not all paramters are required
    

    bool estimateComplexity(int &timeComplexity, int &memComplexity);  // Does nothing right now  
    bool getSolution(double  &PR); // get the solution for PR
    bool getSolution(mex::vector<Factor> &MAR); // get solution for MAR

    bool runInference(int timeLimit); // runs the inference algorithm for a small amount of time. Meant to be called repeatedly
    void printFactors(mex::vector<mex::Factor> flist);


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
    double totalAvailableTime; // the amount of time available for inference
    mex::vector<Factor> flist;
    bool isExact; // true when we've gotten an exact solution
    double dt;
	
    algOptions options;
    
    //
    // Functions
    //
    
    bool initializeDefault(double totalTime);
    bool doLoopyBP();    // Does Loopy BP on factGraph
    bool doGeneralBP();    // Does General BP on factGraph
    bool doIterativeConditioning();
    bool parseCommandOptions(int argc, char** argv);    // Puts command line options in vm
    double  computeVariableOrder(int numTries, double timeLimit);
    bool fitsMBE(const graphModel& gm, const mex::VarOrder& order, const VarSet* cond);
    double solveMBE(const graphModel& gm, const mex::VarOrder& order);
    bool tryExactPR(const graphModel& gm, const mex::VarOrder& order);
    bool gbpPopulateCliques(mex::gbp& _gbp, const mex::VarOrder& order, size_t& ibound, VarSet* cond);
    bool readEvidenceFile();
    mex::vector<Factor> readUaiFile();

  };

}  // namespace lgbp

#endif /* BPINTERFACE_H_ */
