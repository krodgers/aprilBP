/*
 * BPinterface.cpp
 * Belief Propagation interface for Figaro
 *
 *  Created on: 28 Dec 2014
 *      Author: kathryn rodgers
 */

#include "BPinterface.h"
#include "mbe.h"
#include <cstring>
#include <cstdio>
#include <iostream>
#include <fstream>

using namespace std;
using namespace lgbp;
//
// public Functions
//
/* Initialize the computation from options struct. Returns true on success.
// useDefault: if true, uses default option values
// totalTime: total amount of time available to run inference
opts: given options to use; any non specified options should have value -1/NULL
*/
bool BpInterface::initialize(algOptions opts, bool useDefault, double totalTime){
  
  bool result;
  options = algOptions();
  
  // Check files are valid
  ifstream is; is.open(opts.problemFile);
  if (!is.is_open()) {
    if(options.doVerbose)
      std::cout<<"Failed to open problem file"<<std::endl;
    return false;
  }
  is.close();
  if(options.task != Task::PR && options.task != Task::MAR){
    if(options.doVerbose)
      std::cout<<"Missing task"<<std::endl;
    return false;
  }

// set up factors
  if (!readUaiFile()) return false;
  // read evidence
  if(options.evidenceFile != NULL)
    if(!readEvidenceFile())
      return false;
  
 
  // Check all option variables
  // Initializes with a moderate time limit
  options.lbpTime == 0 ?  10 : options.lbpTime; 
  options.lbpIter == 0 ? 2000 : options.lbpIter;
  options.nOrders == 0?  1000 : options.nOrders; 
  options.timeOrder == 0 ?  60: options.timeOrder; 
  options.nExtra == 0 ?  3: options.nExtra;
  options.gbpTime == 0 ?  1000: options.gbpTime; 
  options.gbpIter == 0 ?  -1: options.gbpIter;
  dt == 0 ?  30: dt;
  options.iboundInit == 0 ?  30: options.iboundInit;
  options.doCond == 0 ?  1: options.doCond; 
  options.MemLimit ==0 ? 2*1024.0 : options.MemLimit;

 // set up data structures
  isExact = false;
  phase = Phase::LBP;
  logZ = 0;
  
  
  return true;
}

  
// Initialize the computation using command line arguments
bool BpInterface::initialize(int argc,char** argv )  {
  options = algOptions();
  bool result;
  result = parseCommandOptions(argc, argv);
// set up factors
  if (!readUaiFile()) return false;
  // read evidence
  if(options.evidenceFile != NULL)
    if(!readEvidenceFile())
      return false;

  // set up data structures
  isExact = false;
  phase = Phase::LBP;
  logZ = 0;
  

  return result;
    
}
  
// Initialize using default parameters
bool BpInterface::initialize(double totalTime, char* task, char* problemFile, char* orderFile, char* evidenceFile, bool verbose){
  options = algOptions();
  if ( problemFile == NULL) {
    if(verbose)
      std::cout<< "Missing problem file" << std::endl;
    return false;
  }
  // Check problem file
  ifstream is; is.open(problemFile);
  if (!is.is_open()) {
    if(options.doVerbose)
      std::cout<<"Failed to open problem file"<<std::endl;
    return false;
  }
  is.close();
  // if (strcmp(task, "PR") != 0 && strcmp(task, "MAR") !=0){
  //   if(verbose) std::cout<< "Missing task\n";
  //   return false;
  // }
  options.doVerbose = verbose;
  options.problemFile = problemFile;
  options.orderFile = orderFile;
  options.evidenceFile = evidenceFile;
  options.task = Task(task);

  // set up factors
  if (!readUaiFile()) return false;
  // read evidence
  if(evidenceFile != NULL)
    if(!readEvidenceFile())
      return false;


  if (options.MemLimit <= 0) 
    options.MemLimit = 2048.0;
  mex::randSeed( 1234 );
  if (totalTime <= 0) 
    totalAvailableTime = 30;
  if (totalAvailableTime < 60.0) {								
    // Small time limit
    options.lbpTime = 1.5;  options.lbpIter = 2000;  
    options.nOrders = 100; options.timeOrder = 1.5; options.nExtra = 2;
    options.gbpTime = 300; options.gbpIter = -1; dt = 0.5;
    options.iboundInit = 30; options.MemLimit = std::min(options.MemLimit, 300.0);
    options.doCond = 1;options.ijgp = false;
  }	else if (totalAvailableTime < 1800) {					// Moderate time limit
    options.lbpTime = 10;  options.lbpIter = 2500;
    options.nOrders = 1000; options.timeOrder = 60; options.nExtra = 3;
    options.gbpTime = 1000; options.gbpIter = -1; dt = 30;
    options.iboundInit = 30;
    options.doCond = 1; options.ijgp = true;
  } else {																// Large time limit
    options.lbpTime = 15;  options.lbpIter = 3000;
    options.nOrders = 1000; options.timeOrder = 100; options.nExtra = 3;
    options.gbpTime = 1000;  options.gbpIter = -1; dt = 30;
    options.iboundInit = 30; 
    options.doCond = 1; options.ijgp = true;
  }
  if (options.task==Task::PR) options.gbpIter=1;
  if (options.task==Task::MAR) options.gbpIter=2;

  if(options.doVerbose) std::cout<<"Memory limit set to "<<options.MemLimit<<"mb\n";

  // set up data structures
  isExact = false;
  phase = Phase::LBP;
  logZ = 0;
  flag = Phase::LBP;
  
  return true;
}


// 
bool BpInterface::estimateComplexity(int &timeComplexity, int &memComplexity){
  return true;
}

/*
  Runs inference for the allotted time
*/
bool BpInterface::runInference()
{
  if(isExact || flag == Phase::DONE || phase == Phase::DONE) // answer already exact
    return true;

  //TODO:: check flag before setting phase
  //TODO:: init flag
  double startTime = timeSystem();
  bool success;

  
  switch(phase){
    
  case Phase::LBP: 
    if(flag == Phase::DONE || isExact)
      return true;
    if(!doLoopyBP()) return false;
    phase = Phase::GBP;
    
  case Phase::GBP:
    if(flag == Phase::DONE || isExact)
      return true;
    if(! doGeneralBP()) return false;
    phase = Phase::ITERCOND;
    
  case Phase::ITERCOND:  
    if(flag == Phase::DONE || isExact)
      return true;
    if(!doIterativeConditioning()) return false;
    phase = Phase::LBP;
  } 
  return success;
}
/*
  returns the solution

  Puts the solution in the given argument
  returns true if the solution is successfully given
*/
bool BpInterface::getSolution(mex::vector<Factor> &MAR){
  if(task == Task::PR){
    if(options.doVerbose)
      printf("Failed to get solution -- wrong task");
    return false;
  }
  if(MAR.size() != bel->size()){
    if(options.doVerbose)
      printf("Failed to get solution -- MAR input param wrong size");
    return false;
  }

  mex::vector<Factor>::iterator MARiter;
  mex::vector<Factor>::iterator belIter;

  // ToDO:: FIX THIS --> write out the right thing
  //void writeMAR(const char* outfile, mex::vector<Factor>& fs) {
  // ofstream os(outfile);
  // //os<<"MAR\n1\n";		// !!! 2011 PIC version: need "1" evidence
  // os<<"MAR\n";
  // os<<fs.size()<<" ";
  // for (size_t f=0;f<fs.size();++f) {
  //   os<<fs[f].nrStates()<<" ";
  //   for (size_t i=0;i<fs[f].nrStates();++i) os<<fs[f][i]<<" ";
  // }
  // os<<"\n";
  // os.close();
  // std::cout<<"Wrote MAR\n";


  for(MARiter = MAR.begin(), belIter = bel->begin(); MARiter != MAR.end(); MARiter++, belIter++){
    *MARiter = *belIter;
  }

  return true;
  
}

/*
  returns the solution

  Puts the solution in the given argument
  returns true if the solution is successfully given
*/
bool BpInterface::getSolution(double  &PR)
{
  if(task == Task::MAR){
    if(options.doVerbose)
      printf("Failed to get solution -- wrong task");
    return false;
  }

  PR = logZ;
  return true;
}

/*
  Prints out the factors given in flist
*/
void BpInterface::printFactors(mex::vector<Factor> *flist){
  for (size_t f=0;f<flist->size();++f)  {        
    std::cout << "(*flist)[" <<f<< "]: nvar(" <<(*flist)[f].nvar() << "), numStates("<<  (*flist)[f].nrStates()<< ")"<<std::endl;
    std::cout <<"(*flist)[" << f << "]: values:" ;
    const double *values = (*flist)[f].table();
    for (int v = 0; v < (*flist)[f].nrStates(); v++){
      std::cout<<values[f] << " ";
    }
    printf("\n");
    
  }
}


//
// Private Functions
//

// Parses a command line to initialize options
bool BpInterface::parseCommandOptions(int argc, char** argv){

  po::options_description desc("Available options");
  desc.add_options()
    ("help", "print help message")
    ("file,f", po::value<std::string>(), "input problem filename")
    ("evidence,e", po::value<std::string>(), "input evidence filename")
    ("query,q", po::value<std::string>(), "input marginal map query filename")
    ("seed,S", po::value<int>(),         "random number initial seed")
    ("task,T", po::value<std::string>(), "inference task string")
    ("ibound,i", po::value<int>(&options.iboundInit),       "initial i-bound")
    ("orders,o",    po::value<int>(&options.nOrders)->default_value(1),      "number of variable orderings to try")
    ("order-time,t",po::value<double>(&options.timeOrder)->default_value(1), "max time spend on variable orderings")
    ("order-rand",  po::value<int>(&options.nExtra)->default_value(0),   "var order randomness; n=-1 (none), or among best+n")
    ("order-file", po::value<std::string>(), "problem elimination ordering filename")
    ("memory,m", po::value<double>(&options.MemLimit)->default_value(2*1024.0),    "memory bound (MB)")
    ("lbps", po::value<double>(&options.lbpTime)->default_value(300),  "loopy belief propagation stop (seconds)")
    ("lbpi", po::value<double>(&options.lbpIter)->default_value(2000), "loopy belief propagation stop (iterations)")
    //    ("lbpo", po::value<double>(&options.lbpObj )->default_value(-1),   "loopy belief propagation stop (objective)")
    ("gbps", po::value<double>(&options.gbpTime)->default_value(300),  "gen belief propagation stop (seconds)")
    ("gbpi", po::value<double>(&options.gbpIter)->default_value(-1),   "gen belief propagation stop (iterations)")
    //    ("gbpo", po::value<double>(&options.gbpObj )->default_value(-1),   "gen belief propagation stop (objective)")
    ("condition", po::value<int>(&options.doCond )->default_value(0),"do conditioning search after gbp (0=no,1=incremental,k=at most k states)")
    ("verbose", "verbose output during algorithm execution")
    ("ijgp", "use ijgp regions only")
    ;

  // Add positional defaults for basic UAI14 competition calls
  // po::positional_options_description p;
  // p.add("file", 1);
  // p.add("evidence", 1);
  // p.add("query", 1);
  // p.add("task", 1);

  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  //  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm); 
  po::notify(vm);
  const  char* taskName;
  /*** ARGUMENT CHECKING *************************************************************/
  if (vm.count("help")) { std::cout<<desc<<"\n"; return false; }
  if (vm.count("verbose")) options.doVerbose=true; else options.doVerbose=false;
  if (vm.count("file")) { options.problemFile=vm["file"].as<std::string>().c_str(); }
  else { std::cout<<"Missing input problem file!\n"; return false; }
  if (vm.count("seed")) { mex::randSeed( vm["seed"].as<int>() ); }
  if (vm.count("task")) { taskName= vm["task"].as<std::string>().c_str(); task=Task(taskName); }
  else { std::cout<<"Missing task!\n"; return false; }
  return true;
} 


/*
  reads in an evidence file
*/
bool BpInterface::readEvidenceFile(){
  /*** READ IN EVIDENCE FILE *********************************************************/
  ifstream is2;
  is2.open(options.evidenceFile);
  if (is2.is_open()) {
    std::cout<<"Got evidence file "<<options.evidenceFile<<"\n";
    std::map<uint32_t,size_t> evid;
    int nEvidVar; is2 >> nEvidVar;
    for (size_t i=0;i<nEvidVar;i++) {
      uint32_t vid; size_t vval; is2>>vid>>vval; 
      evid[vid]=vval; evVar |= Var(vid,0);
      //xhat[vid]=vval;
      
      std::cout<<"Evidence on variables "<<evVar<<"\n";
      for (size_t f=0;f<facts->size();f++) {
        if ((*facts)[f].vars().intersects(evVar)) {
          VarSet overlap = (*facts)[f].vars() & evVar;
          evVar += overlap;         // correct missing dimension information
          for (size_t v=0;v<overlap.nvar();++v) 
            (*bel)[overlap[v].label()]=Factor::delta(overlap[v],evid[overlap[v].label()]);
          (*facts)[f] = (*facts)[f].condition( overlap, sub2ind(overlap,evid) );
	  // !!! TODO:FIX: if no variables left, create delta functions with constant value?
        }
      }
      // !!! TODO:DEBUG: if no variables left, create delta functions with constant value?
      //for (size_t i=0;i<nEvidVar;++i) facts.push_back( Factor::delta(evVar[i],evid[evVar[i]]) );
    }
  } else std::cout<<"Evidence file not specified or not found\n";
  
  
}

/*
  Reads in a uai file to facts
  returns true if successful
*/

bool BpInterface::readUaiFile()
{
  /*** READ IN PROBLEM FILE **********************************************************/
  std::cout<<"Reading model file: "<<options.problemFile<<"\n";
  ifstream is; is.open(options.problemFile);
  if(!is.is_open())
    return false;
  
  facts = new mex::vector<Factor>(Factor::readUai10(is));
  if(facts == NULL) return false;
  nvar=0;
  for (size_t f=0;f<facts->size();++f)  {                      // find maximum variable label
    if ((*facts)[f].nvar() > 0){
      nvar=std::max(nvar,(size_t)((*facts)[f].vars().rbegin()->label()+1));
    }
  }

  /////////////// DELETE ME ////////////////
  // printf("\nInitialized Factors: \n");
  // printFactors(facts);
  //printf("\n\n");
  /////////////////////////////////
  bel = new mex::vector<Factor>(nvar);
  return true;
  
}
double  BpInterface::computeVariableOrder(int numTries, double timeLimit){
  // start with a random order in case the model is very dense
  double memUseRandom;
  if(totalAvailableTime < 60.0) memUseRandom = std::exp(28);
  else if (totalAvailableTime < 1800.0)  memUseRandom=std::exp(30);
  else memUseRandom=std::exp(40);

  
  double score = factGraph.order( mex::graphModel::OrderMethod::Random, order, 0, memUseRandom );
  if (options.MemLimit > 0) {
    //    const char *orderFile = NULL;				// Check for pre-specified elimination order file
    ifstream orderIStream; if (options.orderFile!=NULL) orderIStream.open(options.orderFile);
    if (orderIStream.is_open()) { 
      // If we were given an input file with an elimination ordering, just use that
      std::cout << "Reading elimination order from "<<options.orderFile<<"\n";
      size_t ordersize;  orderIStream>>ordersize; 
      assert(ordersize == order.size());
      for (size_t i=0;i<order.size();++i) { size_t tmp;  orderIStream>>tmp; order[i]=tmp;  };
      orderIStream.close();
    } else {
      // Otherwise, calculate elimination order(s) ////////////////////////////////////
      double startOrder = timeSystem();
      size_t iOrder = 0;			
      // Try to build new orders until time or count limit reached ////////////////////
      while (iOrder < options.nOrders && (timeSystem()-startOrder < options.timeOrder)) {
	score = factGraph.order(mex::graphModel::OrderMethod::WtMinFill, order, options.nExtra, score);
	++iOrder;
      }
      double InducedWidth;
      if (score < memUseRandom)  InducedWidth = factGraph.inducedWidth(order);
      std::cout<<"Best order of "<<iOrder<<" has induced width "<<InducedWidth<<", score "<<score<<"\n";

      // If we were given an ordering file name but no file, write our order out 
      ofstream orderOStream; if (options.orderFile!=NULL) orderOStream.open(options.orderFile);
      if (orderOStream.is_open()) {
	std::cout << "Writing elimination order to "<<options.orderFile<<"\n";
	orderOStream<<order.size(); 
	for (size_t i=0;i<order.size();++i) orderOStream<<" "<<order[i]; 
	orderOStream<<"\n";
	orderOStream.close();
      }
    }
  }
} 

// Does Loopy BP on factGraph
bool BpInterface::doLoopyBP() {
  /*** LOOPY BELIEF PROPAGATION ******************************************************/

  //////// DELETE ME ///////
  //  options.lbpTime = 5000000;
  ////////////////////////

  mex::graphModel fg(*facts);

  mex::lbp _lbp(*facts); 
  std::cout<<"Model has "<<nvar<<" variables, "<<_lbp.nFactors()<<" factors"<<std::endl;
  if (options.lbpIter != 0 && options.lbpTime > 0) {
    _lbp.setProperties("Schedule=Priority,Distance=L1");
    _lbp.setStopIter(options.lbpIter); 
    _lbp.setStopMsg(-1);
    // _lbp.setStopObj(options.lbpObj);
    _lbp.setStopObj(-1);
    _lbp.setStopTime(options.lbpTime);
    _lbp.init();
  ///////////////////// DELETE ME ////////////////////////
  printf("Beginning Loopy BP\n");
  // printf("BP factors:\n");
  // mex::vector<Factor> temp = _lbp.beliefs();
  // printFactors(&temp);
  //////////////////////////////////////////
    
    _lbp.run();
    switch (task) {
    case Task::PR:
      logZ = _lbp.logZ() / c_log10; 
      break;
    case Task::MAR: {
      for (size_t v=0;v<nvar;++v) 
	if (!evVar.contains(Var(v,0))) 
	  (*bel)[v]=_lbp.belief( _lbp.localFactor(v) );
    } break;
    }
  
    std::cout<<"LBP "<<logZ<<std::endl;

    _lbp.reparameterize();                                         // Convert loopy bp results to model
    factGraph=graphModel(_lbp.factors());
    /////////////////////// DELETE ME ///////////////////////////
    // printf("Loopy BP factors after reparameterization \n");
    // //printFactors(_lbp.factors());
    // printf("FG Factors\n");
    // printFactors(fg.factors());
    //////////////////////

  }
}



// Does General BP on factGraph
bool BpInterface::doGeneralBP() {
  size_t  ibound = options.iboundInit;
  computeVariableOrder(options.nOrders, options.timeOrder);
  // Try an exact solver and quit if it fits in memory
  if (options.task==Task::PR && tryExactPR(factGraph,order)){
    isExact = true;
    return true;   // TODO: Need to debug without this...
  }

  mex::gbp _gbp = mex::gbp(factGraph.factors());                      // Create a GBP object for later use
  if (options.ijgp) _gbp.setMinimal(true);      // use "ijgp-like" regions? (true = remove all with c=0)
  else              _gbp.setMinimal(false);

  bool doneGBP=false, isExact=false;
  if (options.doVerbose) std::cout<<"\n"<<"Beginning basic GBP...\n";
  while (!doneGBP) { 
    if (options.MemLimit > 0) { 
      _gbp.clearRegions();
      isExact = gbpPopulateCliques(_gbp,order,ibound,NULL);
    } else {
      // Didn't build MBE => no way to check ibound
      ibound = 0; isExact=false;
    }

    try {
      // Run GBP on the region graph 
      std::cout<<"GBP with "<<_gbp.nRegions()<<" regions; mem "<<_gbp.memory()<<"M\n";
      _gbp.setProperties("Schedule=Fixed");
      _gbp.init();
      _gbp.setStopIter(options.gbpIter); 
      //_gbp.setStopObj(options.gbpObj); _gbp.setStopMsg(-1.0); 
      _gbp.setStopObj(-1); _gbp.setStopMsg(-1.0); 
      _gbp.setVerbose(options.doVerbose);
      if (isExact) _gbp.setDamping(-1.0);  // no damping if we think it's exact
      // Get region indices for single-variable beliefs
      mex::vector<mex::gbp::findex> regions(nvar);
      for (size_t v=0;v<nvar;++v) {
	if (!evVar.contains(Var(v,0))) regions[v]=_gbp.regionWith(Var(v,0));
      }

      double gbpLeft, gbpStop = timeSystem()+options.gbpTime;
      while ( (gbpLeft = gbpStop - timeSystem()) > 0 ) {
	_gbp.setStopTime( std::min( dt , gbpLeft ) );
	_gbp.run();
	switch (task) {
	case Task::PR: logZ = _gbp.logZStable()/c_log10; break; 
	case Task::MAR: {
	  //for (size_t v=0;v<nvar;++v) if (!evVar.contains(Var(v,0))) bel->at(v)=_gbp.belief(Var(v,0));
	  for (size_t v=0;v<nvar;++v) 
	    if (!evVar.contains(Var(v,0))) 
	      bel->at(v) =_gbp.computeRegionBelief(regions[v]).marginal(Var(v,0));
	} break;
	}
	std::cout<<"GBP "<<_gbp.logZ()/c_log10<<"\n";
	//	if (_gbp.dObj() < options.gbpObj) { std::cout<<"Reached objective tolerance\n"; break; }
	if (_gbp.iter() >= options.gbpIter && options.gbpIter > 0) { std::cout<<"Reached iteration limit\n"; break; }
	if (_gbp.logZ() == -mex::infty()) { std::cout<<"Model deemed inconsistent\n"; break; }
		
      }
      doneGBP = true;

      if (isExact && _gbp.dObj()<.0001) { //InducedWidth <= ibound)
	std::cout<<"Answer should be exact\n";
	return 0;
      }

    } catch (std::exception& e) {
      if (_gbp.getDamping() > 0.25) {
	doneGBP=true; std::cout<<"Caught exception (memory overreach?).  Damping on => quitting GBP\n";
      } else {
    	doneGBP=false; options.MemLimit*=.9; ibound--;
	std::cout<<"Caught exception (memory overreach?).  Trying again with ibound "<<ibound<<" and options.MemLimit "<<options.MemLimit<<"\n";
      }
    }
  }
}



bool BpInterface::doIterativeConditioning(){
  /*** ITERATIVE CONDITIONING AND GBP *************************************/
  // More general way to do this: condition-min (=0), condition-max (=inf)   (# states?)
  // flag to decide whether to discard if non-convergent (?) (see this vs cgbp) (or ok with damping?)
  // decide how much time?  gbpTime, or as a fraction of # states?
  // =1 => previous code; > 1 is subsequent code
  // need to include exact lnZ check here also


  ///////////// DELETE ME /////////////////////////////////
  printf("Iterative Conditioning and GBP\n");

  //////////////////////////////////////////
  bool exact = false;
  size_t ibound = options.iboundInit, InducedWidth=10000;


  if (options.doVerbose) std::cout<<"\n"<<"Beginning conditioned GBP...\n";
  if (order.size()==0) order=factGraph.order(mex::graphModel::OrderMethod::MinWidth);  // need an order if none yet...
 
  // mex::gbp _gbp = mex::gbp( mex::vector<Factor>() );  // !!! blank out GBP object; restore memory
  VarSet cond;

  //!!! TODO: if options.doCond > 1, we should immediately populate with that many states

  while (!isExact) {
    // add check for best conditioner given cardinality limit !!!  or, only increment if anytime...
    cond += factGraph.bestConditioner(order,cond);
    if (options.doVerbose) std::cout<<"\n";
    std::cout<<"Conditioning "<<cond<<"\n";
    ibound = options.iboundInit;		// check all iBounds again (in case higher available)
    bool doneCGBP=false;
    while (!doneCGBP) {
      bool useMBE=false;
      if (task == Task::PR && fitsMBE(factGraph,order,&cond)) {
        std::cout<<"Trying exact via MBE\n";
        useMBE=true; isExact=true;
      }
      try {
    	Factor lnZ(cond);
	bool failed = false;
    	mex::vector<mex::vector<Factor> > condMarginals;
     	mex::vector<mex::gbp::findex> regions(nvar);			// !!! dangerous; pull out & match to prev (clear+resize)
    	for (size_t i=0;i<lnZ.nrStates();++i) {
	  std::map<Var,size_t> val;  ind2sub(cond,i,val);
	  mex::vector<Factor> fcond = factGraph.factors();
	  for (size_t f=0;f<fcond.size();++f) {
	    VarSet isect = cond & fcond[f].vars();
	    if (isect.size() > 0) fcond[f] = fcond[f].condition(isect, sub2ind(isect,val));
	  }

	  try { if (useMBE) {         // if a bucket elim pass was good enough, do that:
	      mex::graphModel gm(fcond);
	      lnZ[i] = solveMBE(gm,order);
	      for (size_t v=0;v<cond.size();++v) std::cout<<cond[v]<<"="<<val[cond[v]]<<" "; std::cout<<lnZ[i]<<"\n";
	      continue;
	    }                     // otherwise we need to do GBP-like updates:
	  } catch (std::exception& e) {
	    std::cout<<"Caught exception; failure in MBE; trying GBP\n";
	  }

	  mex::lbp fgcond(fcond); fgcond.setProperties("Schedule=Priority,Distance=L1");
	  fgcond.setStopIter(options.lbpIter); fgcond.setStopMsg(-1); 
	  //  fgcond.setStopObj(options.lbpObj);
	  fgcond.setStopObj(-1);
	  fgcond.setStopTime(0.5);
	  fgcond.init();
	  fgcond.run();
	  fgcond.reparameterize();
	  fcond = fgcond.factors();

	  // TODO:DEBUG: it seems that reinitializing gbp is not equivalent to re-constructing a new one.. FIX
	  mex::gbp _gbp(fcond);
	  if (options.ijgp) _gbp.setMinimal(true); else _gbp.setMinimal(false);  // use "ijgp-like" regions?
	  isExact = gbpPopulateCliques(_gbp,order,ibound,&cond);
	  _gbp.setProperties("Schedule=Fixed");
	  _gbp.setStopIter(options.gbpIter); 
	   _gbp.setStopMsg(-1.0);
	   //_gbp.setStopObj(options.gbpObj);
	   _gbp.setStopObj(-1);
	  _gbp.setStopTime(options.gbpTime); _gbp.setVerbose(options.doVerbose);
	  _gbp.init();

	  if (task==Task::MAR) {		// TODO: change to loop over "inferred" variable list...
	    for (size_t v=0;v<nvar;++v) {
	      if (!evVar.contains(Var(v,0)) && !cond.contains(Var(v,0))) regions[v]=_gbp.regionWith(Var(v,0));
	    }
	  }
	  // uncomment duplicate lines below

	  // !!! make incremental on last output round
	  //      _gbp.setFactors(fcond); _gbp.setProperties("Schedule=Fixed");
	  //      _gbp.init();

	  _gbp.run();
	  //if (_gbp.dObj() >= gbpObj) { failed=true; break; }     // convergence failure on this condition
	  //TODO: never worry about this?  or, only worry if anytime enabled?
	  lnZ[i] = _gbp.logZStable();
	  if (task==Task::MAR) {
	    condMarginals.push_back(*bel);
	    for (size_t v=0;v<nvar;++v)
	      if (!evVar.contains(Var(v,0)) && !cond.contains(Var(v,0)))
            	condMarginals[i][v]=_gbp.computeRegionBelief(regions[v]).marginal(Var(v,0));
	    //condMarginals[i][v]=_gbp.belief(Var(v,0));
	  }
	  for (size_t v=0;v<cond.size();++v) std::cout<<cond[v]<<"="<<val[cond[v]]<<" "; std::cout<<lnZ[i]<<"\n";
    	}
    	if (failed) {std::cout<<"Failing out\n"; doneCGBP=true; continue;} // if we failed out, condition on more vars
    	
	doneCGBP=true;
    	double lnZtot = lnZ.logsumexp();
    	switch (task) {
      	case Task::PR:
	  logZ = lnZtot/c_log10;
	  break;
      	case Task::MAR:
	  Factor probs = (lnZ - lnZtot).exp();
	  for (size_t v=0;v<nvar;++v) {
	    if (evVar.contains(Var(v,0))) { } // evidence variables not updated
	    else if (cond.contains(Var(v,0)))  { bel->at(v) = probs.marginal(Var(v,0)); }
	    else {		// TODO: change to incremental update?
	      bel->at(v) = condMarginals[0][v] * probs[0];
	      for (size_t i=1;i<lnZ.nrStates();++i) bel->at(v) += condMarginals[i][v] * probs[i];
	    }
	  }
	  break;
    	}
    	std::cout<<"Conditioning "<<cond<<" => "<<lnZtot<<" ("<<lnZtot/c_log10<<")\n";

      } catch (std::exception& e) {
      	doneCGBP=false; options.MemLimit*=.9; ibound--;
      	std::cout<<"Caught exception (memory overreach?).  Trying again with ibound "<<ibound<<" and options.MemLimit "<<options.MemLimit<<"\n";
      	continue;
	// TODO: this is right if we're only doing this part, but not right if we're running incrementally
      }

      // !!! TODO: if options.doCond > 1, quit (non-incremental)?
      if (isExact) { std::cout<<"Answer should be exact\n"; return 0; }
    }
  }

  // run CCS expansion on larger regions until timeout (?)
}


//////////////////////////////////
// Helper Functions 
//////////////////////////////

bool BpInterface::fitsMBE(const graphModel& gm, const mex::VarOrder& order, const VarSet* cond) {
  double mbCutoff = options.MemLimit/sizeof(double)*1024*1024;     // translate memory into MBE cutoff
  bool isExact = true;
  mex::mbe mb(gm.factors());
  mb.setOrder(order);
  mb.setProperties("ElimOp=SumUpper,sBound=inf,DoMatch=1,DoMplp=0,DoFill=0,DoJG=0,DoHeur=0");
  mb.setIBound(100); double mbMem = mb.simulateMemory(NULL,cond,mbCutoff,&isExact);
  if (mbMem < mbCutoff && isExact) return true;
  return false;
}

double BpInterface::solveMBE(const graphModel& gm, const mex::VarOrder& order) {
  mex::mbe mb(gm.factors());
  mb.setOrder(order);
  mex::vector<mex::mbe::vindex> pt(order.size()); for (size_t v=0;v<pt.size();++v) pt[v]=-1;
  mb.setPseudotree(pt);
  mb.setProperties("ElimOp=SumUpper,sBound=inf,DoMatch=1,DoMplp=0,DoFill=0,DoJG=0,DoHeur=0");
  mb.setIBound(100); //double mbMem = mb.simulateMemory(NULL,NULL,mbCutoff,&isExact);
  std::cout<<"Attempting exact solve\n";
  //std::cout<<"Attempting exact solve; mbMem="<<mbMem<<" vs "<<mbCutoff<<" ("<<options.MemLimit<<")\n";
  mb.init();
  return mb.logZ();
}


bool BpInterface::tryExactPR(const graphModel& gm, const mex::VarOrder& order) {

  ///////////////// DELTE ME ////////////////////////////
  printf("tryExactPR:\n");

  ////////////////////////
  try {
    double mbCutoff = options.MemLimit/sizeof(double)*1024*1024;     // translate memory into MBE cutoff
    bool isExact = true;
    mex::mbe mb(gm.factors());
    
    ////////////////DELETE ME////////////////////
    //order[0] = 0; order[1] = 1;
    // printf("MBE Factors\n");
    // mex::vector<Factor> temp = mb._gmo.factors();
    // printFactors(&temp);
    // std::cout << "Order: "<< order[0] << " " <<order[1]<< std::endl;
    /////////////////////////////////////
    
    mb.setOrder(order);
    mb.setProperties("ElimOp=SumUpper,sBound=inf,DoMatch=1,DoMplp=0,DoFill=0,DoJG=0,DoHeur=0");
    mb.setIBound(100); double mbMem = mb.simulateMemory(NULL,NULL,mbCutoff,&isExact);
    if (mbMem < mbCutoff && isExact) {
      std::cout<<"Attempting exact solve; mbMem="<<mbMem<<" vs "<<mbCutoff<<" ("<<options.MemLimit<<")\n";
      mb.init();
      logZ = mb.logZ()/c_log10;
      std::cout<<"Exact solution by MBE: "<<logZ<<"\n";
      return true;
    }
    
  } catch (std::exception& e) {
    // Failed (probably for memory reasons) => try GBP
    std::cout<<"Failed (due to memory problem?)  Trying GBP\n";
  }
  return false;
}


bool BpInterface::gbpPopulateCliques(mex::gbp& _gbp, const mex::VarOrder& order, size_t& ibound, VarSet* cond) {
  ///////////////////// DELETE ME ////////////////////
  printf("gbpPopulateCliques\n");
  ///////////////////////////////////////
  bool isExact = false;
  double mbCutoff = options.MemLimit/sizeof(double)*1024*1024;     // translate memory into MBE cutoff
  _gbp.clearRegions();
  mex::mbe mb(_gbp.gmOrig().factors());
  mb.setOrder(order);
  mex::vector<mex::mbe::vindex> pt(order.size()); for (size_t v=0;v<pt.size();++v) pt[v]=-1;
  mb.setPseudotree(pt);
  mb.setProperties("ElimOp=SumUpper,sBound=inf,DoMatch=1,DoMplp=0,DoFill=0,DoJG=1,DoHeur=0");
  mb.setIBound(ibound);
  mex::vector<VarSet> cliques;
  double mem = std::numeric_limits<double>::infinity();
  mbCutoff *= 2;  // leave a little slack for mismatch in criteria
  double mbMem = mb.simulateMemory(&cliques,cond, mbCutoff, &isExact);
  if (mbMem < mbCutoff) { _gbp.addRegions(cliques); mem = _gbp.memory(); }
  while ( ibound > 0 && (mbMem >= mbCutoff || mem > options.MemLimit) ) {
    std::cout<<"MBE iBound "<<ibound<<" = "<<mem<<"M\n";
    mb.setIBound(--ibound); cliques.clear(); mbMem=mb.simulateMemory(&cliques,cond, mbCutoff, &isExact);
    if (mbMem < mbCutoff) { _gbp.clearRegions(); _gbp.addRegions(cliques); mem=_gbp.memory(); }
  }
  //ofstream ofs("cliques.mbe.txt");
  //for (size_t c=0;c<cliques.size();++c) ofs<<cliques[c]<<"\n"; std::cout<<"\n";  // output for DEBUG !!!
  //ofs.close();
  std::cout<<"MBE iBound "<<ibound<<" = "<<mem<<"M"<<std::endl;
  return isExact;
}

