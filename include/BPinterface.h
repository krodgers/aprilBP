#ifndef BPINTERFACE_H_
#define BPINTERFACE_H_
 
#include <stdlib.h>
#include <stdio.h>
#include <sstream>

#include <string>
#include <list>

//#include "boost/program_options.hpp"

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

//namespace po = boost::program_options;


namespace lgbp {
#define c_log10 std::log(10)
  typedef struct {
    double MemLimit;
    double lbpTime, lbpIter;
    double gbpTime, gbpIter;
    int    doCond;
    bool   doVerbose; // toggles printing to std::out
    int    iboundInit; // initial ibound liimit
    double timeOrder; // max time to spend on variable ordering
    int    nOrders; // max number of variable orderings to try
    const char* problemFile; // the name of the file to read the problem from
    int task; // which task we're doing
    const char*  orderFile; // name of the file to read a variable ordering from
    const char* evidenceFile; // name of the file to read evidence from
    // std::string query; // name of the file to read MMAP query from
    int  nExtra;
    bool ijgp;
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
    BpInterface();
	
    ~BpInterface();
    
    //   bool initialize(int argc, char** argv); // Initialize the computation. Returns true on success.
    bool initialize(algOptions opts, bool useDefault, double totalTime); // Initialize with given paramters; not all paramters are required
    bool initialize(double totalTime, char* task, char* problemFile, char* orderFile, char* evidenceFile, bool verbose);
 

    bool estimateComplexity(int &timeComplexity, int &memComplexity);  // Does nothing right now  
    bool getSolution(double  &PR); // get the solution for PR
    bool getSolution(mex::vector<Factor> &MAR); // get solution for MAR

    bool runInference(); // runs the inference algorithm
    bool resetSolver();
    bool stopComputation();

  private:
    //
    // Fields
    //
    MEX_ENUM(Phase ,  LBP, GBP, ITERCOND, DONE );

    int phase; // track which part of the inference we're doing
    int flag; // track if inference has been told to stop
    mex::vector<Factor>* bel; // stores results for MAR
    mex::vector<Factor>* facts; // stores the factors 
    VarSet evVar; // evidence variables
    mex::VarOrder order; // the ordering of the variables
    size_t nvar;
    mex::graphModel factGraph; // the graph of the model
    double logZ; // result of PR
    double totalAvailableTime; // the amount of time available for inference
    bool isExact; // true when we've gotten an exact solution
    double dt;
    std::string  logFileName;

    algOptions options;
#ifdef LOGFILE
    std::ofstream *out;
    std::streambuf *restoreSTDCOUT;
#endif
    
    
    //
    // Functions
    //
    bool doLoopyBP();    // Does Loopy BP on factGraph
    bool doGeneralBP();    // Does General BP on factGraph
    bool doIterativeConditioning();
    // bool parseCommandOptions(int argc, char** argv);    // Puts command line options in vm
    double  computeVariableOrder(int numTries, double timeLimit);
    bool fitsMBE(const graphModel& gm, const mex::VarOrder& order, const VarSet* cond);
    double solveMBE(const graphModel& gm, const mex::VarOrder& order);
    bool tryExactPR(const graphModel& gm, const mex::VarOrder& order);
    bool gbpPopulateCliques(mex::gbp& _gbp, const mex::VarOrder& order, size_t& ibound, VarSet* cond);
    bool readEvidenceFile();
    bool readUaiFile();
    void printFactors(mex::vector<mex::Factor> *flist);
    void writeLog(std::string logMsg); // writes logMsg to logfile
    void writeLog(std::stringstream logMsg);
    void writePR(const char* outfile, double logZ);
    void writeMAR(const char* outfile, mex::vector<Factor>& fs);

  };
}  // namespace lgbp

#endif /* BPINTERFACE_H_ */
