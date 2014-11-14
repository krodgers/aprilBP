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

//using namespace mex;
using mex::mxObject;
using mex::Var;
using mex::VarSet;
using mex::Factor;
using mex::vector;
using mex::graphModel;
using mex::factorgraph;
using mex::mplp;
//
using mex::timeSystem;


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
	BpInterface();


  // Initialize the computation. Returns true on success.
	bool BpInterface::initialize(int argc, char** argv, int theTask);

// Does nothing right now  
  bool estimateComplexity();

  // 
  double solvePR();

  double solveMAR();

  double getPRSolution();
  mex::vector<Factor> getMARSolution();

private:
	//
	// Fields
	//
	int task; // which task we're doing
	mex::vector<Factor> bel;
	mex::vector<uint32_t>	xhat;
	mex::vector<Factor> facts;
	mex::graphModel factGraph;
	
	VarSet maxVars;
	
	double MemLimit;
	double lbpTime, lbpIter, lbpObj, lbpErr;
	double gbpTime, gbpIter, gbpObj;
	double dt;
	int    doCond;
	bool   doVerbose; // toggles printing to std::out
	int    nExtra;
	int    iboundInit;
	double logZ;
	
	double timeOrder;
	int    nOrders;
	double memUseRandom = std::exp(40.0);	// if memory *way* out of reach, just use random ordering...
	
	po::variables_map vm; // algorithm options
	
	//
	// Functions
	//
	bool doLoopyBP();
	bool doGeneralBP();
	bool doIterativeConditioning();
	bool parseCommandOptions(int argc, char** argv);
	mex::VarOrder computeVariableOrder(int numTries, double timeLimit);
	bool fitsMBE(const graphModel& gm, const mex::VarOrder& order, const VarSet* cond);
	double solveMBE(const graphModel& gm, const mex::VarOrder& order);
	bool tryExactPR(const graphModel& gm, const mex::VarOrder& order);
	bool gbpPopulateCliques(mex::gbp& _gbp, const mex::VarOrder& order, size_t& ibound, VarSet* cond);
	bool readEvidenceFile(mex::vector<Factor> &flist);
	mex::vector<Factor> readUaiFile();


			
	

  };

}  // namespace daoopt

#endif /* BPINTERFACE_H_ */
