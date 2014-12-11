/*
  Belief propagation interface for APPRIL


  author: Kathryn Rodgers
*/
#include "BPinterface.h"

#include "stdio.h"
#include <iostream>
#include <fstream>

#include "mbe.h"

#define c_log10 std::log(10)

using namespace lgbp;
using namespace mex;
using namespace std;

void BpInterface::printFactors(mex::vector<Factor> flist){
  for(int f = 0; f < flist.size(); f++){
    printf("flist[%d]: nvar(%d), numStates(%d)\n", f, flist[f].nvar(), flist[f].nrStates());
    printf("flist[%d]: values:", f) ;
    const double *values = flist[f].table();
    for (int v = 0; v < flist[f].nrStates(); v++){
      printf("%f  ", values[f]);
    }
    printf("\n");
  }
}


/*
  Initializes parameters and data structures

  Returns true if successful and false otherwise
 */
bool BpInterface::initialize(int argc, char** argv){
  // set up algorithm options
  bool allGood = parseCommandOptions(argc, argv);
  if(!allGood)
    return false;

  // read  in uai file
  facts = readUaiFile();
  if(vm.empty()){
    if(doVerbose)
      printf("Could not read UAI file\n");
    return false;
  }
  allGood = readEvidenceFile();
  if(!allGood && vm.count("evidence"))
    return false;
  
  // set up data structures
  factGraph = mex::graphModel(facts);

  // initialize other member variables
  memUseRandom = std::exp(40.0);
  isExact = false;
  phase = Phase::LBP;
  logZ = 0;
  bel = mex::vector<Factor>();
  
  //// DELETE ME
  printf("Factors:\n");
  printFactors(facts);
  lbpTime = 4;
  gbpTime = 0;
  nOrders = 0;
  doCond = 0;
  ////////////////////////////
  return true;

}

/*
  Read in UAI file
*/
mex::vector<Factor> BpInterface::readUaiFile(){
  const char* probName;
	
  if (vm.count("file")) {
    probName = vm["file"].as<std::string>().c_str();
  }
  else {
    std::cout << "Missing input problem file!\n";
    return mex::vector<Factor>(); // return empty vector
  }
	
  //std::cout << "Reading model file: " << probName << "\n";
  ifstream is; 
  is.open(probName);
  if (!is.is_open())
    throw std::runtime_error("Failed to open problem file");
  mex::vector<Factor> flist = Factor::readUai10(is);
  nvar = 0;
  for (size_t f = 0; f < flist.size(); ++f){		// find maximum variable label
    if (flist[f].nvar() > 0)
      nvar = std::max(nvar, (size_t)(flist[f].vars().rbegin()->label() + 1));
  }
  bel.resize(nvar);
  return flist;
}

/*
  Read in Evidence File
*/
bool BpInterface::readEvidenceFile(){
  ifstream is2;
  if (vm.count("evidence")) {
    is2.open(vm["evidence"].as<std::string>().c_str());
    // Check for empty file
    int c = is2.peek();  // peek character
    if ( c == EOF )
      return false;
  }
  if (is2.is_open()) {
    if (doVerbose)
      std::cout << "Got evidence file " << vm["evidence"].as<std::string>().c_str() << "\n";
    std::map<uint32_t, size_t> evid;
    int nEvid = 1;			// 2014 format: single evidence, no # of evidences entry
    if (nEvid > 0) {
      int nEvidVar; is2 >> nEvidVar;
      for (size_t i = 0; i < nEvidVar; i++) {
	uint32_t vid; size_t vval; is2 >> vid >> vval;
	evid[vid] = vval; evVar |= Var(vid, 0);
      }
      if (doVerbose)
	std::cout << "Evidence on variables " << evVar << "\n";
      for (size_t f = 0; f < facts.size(); f++) {
	if (facts[f].vars().intersects(evVar)) {
	  VarSet overlap = facts[f].vars() & evVar;
	  evVar += overlap;         // correct missing dimension information
	  for (size_t v = 0; v < overlap.nvar(); ++v)
	    bel[overlap[v].label()] = Factor::delta(overlap[v], evid[overlap[v].label()]);
	  facts[f] = facts[f].condition(overlap, sub2ind(overlap, evid));
	  // !!! TODO:FIX: if no variables left, create delta functions with constant value?
	}
      }
      // !!! TODO:DEBUG: if no variables left, create delta functions with constant value?
      //for (size_t i=0;i<nEvidVar;++i) facts.push_back( Factor::delta(evVar[i],evid[evVar[i]]) );
    }
    return true;
  }
  else {
    std::cout << "Evidence file not specified or not found\n";
    return false;
  }

}

/*
  Parses and initialized command-line-given arguments
  Non-positional arguments must be given via --argument, i.e --help
  positional arguments: -f problem_file -e evidence_file -q query-T task 


  returns true if all required arguments are sent and all arguments are successfully parsed
  false otherwise
*/
bool BpInterface::parseCommandOptions(int argc, char** argv){
  // parse command options
  po::options_description desc("Available options");
  desc.add_options()
    ("help", "print help message")
    ("file,f", po::value<std::string>(), "input problem filename")
    ("evidence,e", po::value<std::string>(), "input evidence filename")
    ("query,q", po::value<std::string>(), "input marginal map query filename")
    ("seed,S", po::value<int>(), "random number initial seed")
    ("task,T", po::value<std::string>(), "inference task string")
    ("ibound,i", po::value<int>(&iboundInit), "initial i-bound")
    ("orders,o", po::value<int>(&nOrders)->default_value(1), "number of variable orderings to try")
    ("order-time,t", po::value<double>(&timeOrder)->default_value(1), "max time spend on variable orderings")
    ("order-rand", po::value<int>(&nExtra)->default_value(0), "var order randomness; n=-1 (none), or among best+n")
    ("order-file", po::value<std::string>(), "problem elimination ordering filename")
    ("memory,m", po::value<double>(&MemLimit)->default_value(2 * 1024.0), "memory bound (MB)")
    ("dt", po::value<double>(&dt)->default_value(300), "file write time interval")
    ("lbps", po::value<double>(&lbpTime)->default_value(300), "loopy belief propagation stop (seconds)")
    ("lbpi", po::value<double>(&lbpIter)->default_value(2000), "loopy belief propagation stop (iterations)")
    ("lbpe", po::value<double>(&lbpErr)->default_value(-1), "loopy belief propagation stop (msg error)")
    ("lbpo", po::value<double>(&lbpObj)->default_value(-1), "loopy belief propagation stop (objective)")
    ("gbps", po::value<double>(&gbpTime)->default_value(300), "gen belief propagation stop (seconds)")
    ("gbpi", po::value<double>(&gbpIter)->default_value(-1), "gen belief propagation stop (iterations)")
    ("gbpo", po::value<double>(&gbpObj)->default_value(-1), "gen belief propagation stop (objective)")
    ("condition", po::value<int>(&doCond)->default_value(0), "do conditioning search after gbp (0=no,1=incremental,k=at most k states)")
    ("verbose,v", "verbose output during algorithm execution")
    ("ijgp", "use ijgp regions only")
    ;

  // Add positional defaults for basic UAI14 competition calls
  po::positional_options_description p;

  p.add("file", 1);
  p.add("evidence", 1);
  p.add("query", 1);
  p.add("task", 1);

  // Parse through arguments 	
  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  po::notify(vm);



  /*** ARGUMENT CHECKING *************************************************************/
  if (vm.count("help")) { 
    std::cout << desc << "\n"; 
    return false; }
  if (vm.count("verbose")) 
    doVerbose = true; 
  else 
    doVerbose = false;


  if (!vm.count("file")) { 
    std::cout << "Missing input problem file!\n"; 
    return false;
  }
  if (vm.count("seed")) { 
    mex::randSeed(vm["seed"].as<int>());
  }
  if (vm.count("task")) {
    std::string taskName = vm["task"].as<std::string>();
    task = Task(taskName.c_str());
  }
  else { 
    std::cout << "Missing task!\n"; 
    return false;
  }

  double TotalTime = 60.0;
  const char *inftime = ::getenv("INF_TIME");
  const char *infmem = ::getenv("INF_MEMORY");
  if (infmem != NULL && inftime != NULL) {
    std::cout << "Using default settings given UAI14 environment variables\n";
    if (vm["memory"].defaulted() && infmem != NULL) 
      MemLimit = atof(infmem) * 1000 - 100;
    if (!vm.count("seed")) 
      mex::randSeed(1234);
    if (inftime != NULL) 
      TotalTime = atof(inftime);
    if (TotalTime < 60.0) {	// Small time limit
      lbpTime = 1.5; lbpErr = 1e-4;
      nOrders = 100; timeOrder = 1.5; nExtra = 2;
      gbpTime = 300; gbpObj = 1e-2; gbpIter = -1; dt = 0.5;
      iboundInit = 30; MemLimit = std::min(MemLimit, 300.0); memUseRandom = std::exp(28);
      doCond = 1;
    }
    else if (TotalTime < 1800) {// Moderate time limit
      lbpTime = 10; lbpErr = 1e-4;
      nOrders = 1000; timeOrder = 60; nExtra = 3;
      gbpTime = 1000; gbpObj = 1e-2; gbpIter = -1; dt = 30;
      iboundInit = 30; memUseRandom = std::exp(30);
      doCond = 1;
    }
    else {																// Large time limit
      lbpTime = 15; lbpErr = 1e-4;
      nOrders = 1000; timeOrder = 100; nExtra = 3;
      gbpTime = 1000; gbpObj = 1e-2; gbpIter = -1; dt = 30;
      iboundInit = 30; memUseRandom = std::exp(40);
      doCond = 1;
    }
    if (task == Task::PR) 
      gbpIter = 1;
    if (task == Task::MAR) 
      gbpIter = 2;
  }

  if(doVerbose)
    std::cout << "Memory limit set to " << MemLimit << "mb\n";

  return true;
}


//TODO:: what should this return?  String-Value pairs?
/**
   Gives a shot in the dark guesstimate about the complexity

   parameters:  timeComplexity: out param; where the time complexity is stored (secs)
                memComplexity: out param; where the memory complexity is stored
 */
bool BpInterface::estimateComplexity(int &timeComplexity, int &memComplexity){
  //TODO:: come up with something to return
  // compute induced width of a few orderings and return textbook estimates?	
  double score;
  score = computeVariableOrder(nOrders, 3);
  
  timeComplexity = score;
  memComplexity = factGraph.inducedWidth(order) * nvar;
  
  return true;
}

/**
   Runs the solvers for the specified amount of time.  Can be called repeatedly in succession and the 
   computation will pick up where it left off

   Does NOT return a solution.  Call getXSolution for the solution for task X

   returns true if nothing's gone wrong
 */
bool BpInterface::runInference(int timeLimit){
  // todo:: implement timing

  if(isExact) // already done all the inference possible
    return true;
  bool result;
  // Try exact solution
  //TODO:: where to first compute order?
  switch(phase){
  case Phase::DONE: 
    return true;
  case Phase::LBP:
    // doLoopyBP
    result = doLoopyBP();//TODO:: run more than once? split up time so can stop in the middle?
    phase = Phase::EXACT;
  
  case Phase::EXACT:  
    computeVariableOrder(nOrders, timeOrder);
    
    if (tryExactPR(factGraph, order)){
      //found a solution
      phase = Phase::GBP;
      return getPRSolution();
    }
    break;
  case Phase::GBP:
    // do Gen BP
    if (result) // if there's no problems encountered
      result = doGeneralBP();		
    phase = Phase::ITERCOND;
  case Phase::ITERCOND:
    // do iterCond
    if (result && !isExact)
      result = doIterativeConditioning();
   
    
  }



  
  	
}

/**
   Gets the current solution for PR task

   returns -1 if the task being solved wasn't PR
 */
double BpInterface::getPRSolution(){
  // PR Solution
  if (task != Task::PR){
    if(doVerbose)
      printf("Wrong solution function for task\n");
    return -1;
  }
    
  if (doVerbose)
    std::cout << "Got PR Solution: " << logZ / c_log10 << "\n";
  return  logZ / c_log10;
}

/**
   Gets the current solution for MAR task

   returns empty vector if the task being solved wasn't MAR
 */

mex::vector<Factor> BpInterface::getMARSolution(){
  // MAR solution

 if (task != Task::MAR){
    if(doVerbose)
      printf("Wrong solution function for task\n");
    return mex::vector<Factor>();
  }
 
  // fs: mex::vector<Factor>& fs
  // os << "MAR\n";
  // os << fs.size() << " ";
  // for (size_t f = 0; f < fs.size(); ++f) {
  // 	os << fs[f].nrStates() << " ";
  // 	for (size_t i = 0; i < fs[f].nrStates(); ++i)
  // 		os << fs[f][i] << " ";
  // }
  // os << "\n";
  // os.close();
  // if (doVerbose)
  // 	std::cout << "Wrote MAR\n";
  
  // TODO:: takes a query?

  return bel; 
	
}

/*
  Does Loopy BP 
  Puts results into BpInterface.bel
  Initializes factGraph for use later

  returns if it sucessfully completes
*/
bool BpInterface::doLoopyBP(){
  /*** LOOPY BELIEF PROPAGATION ******************************************************/
  mex::lbp _lbp(facts);

  if (doVerbose)
    std::cout << "Model has " << nvar << " variables, " << _lbp.nFactors() << " factors\n";
  if (lbpIter != 0 && lbpTime > 0) { // TODO:: change if to while?
    _lbp.setProperties("Schedule=Priority,Distance=L1");
    _lbp.setStopIter(lbpIter);
    _lbp.setStopMsg(lbpErr);
    _lbp.setStopObj(lbpObj);
    _lbp.setStopTime(lbpTime);
    _lbp.init();
		
    _lbp.run();
  
    switch (task) {
    case Task::PR: 
      logZ = _lbp.logZ();
      break;
    case Task::MAR: {
      for (size_t v = 0; v < nvar; ++v) {
	if (!evVar.contains(Var(v, 0)))
	  bel[v] = _lbp.belief(_lbp.localFactor(v));
      }
    }break;
    }
  
    if (doVerbose)
      std::cout << "LBP " << _lbp.logZ() / c_log10 << "\n";
    _lbp.reparameterize();                                         // Convert loopy bp results to model
    factGraph = graphModel(_lbp.factors());
    printFactors(_lbp.factors());
    
  }		
  
  return true;
   
}
/*
  Runs General Belief Propagation 

  Puts results in bel or logZ
  returns true if successfully completes
*/
bool BpInterface::doGeneralBP() {
  bool res = true;
  mex::gbp _gbp(factGraph.factors());                      // Create a GBP object for later use
  size_t ibound = iboundInit, InducedWidth=10000;
  double score;

  score = computeVariableOrder(nOrders,timeOrder);
  if (vm.count("ijgp"))
    _gbp.setMinimal(true);      // use "ijgp-like" regions? (true = remove all with c=0)
  else
    _gbp.setMinimal(false);
  
  bool doneGBP = false, isExact = false;
  while (!doneGBP) {
    if (MemLimit > 0) {
      _gbp.clearRegions();
      isExact = gbpPopulateCliques(_gbp, order, ibound, NULL);
    }
    else {
      // Didn't build MBE => no way to check ibound
      ibound = 0; isExact = false;
    }
    try {
      // Run Gbp On The Region Graph
      if (doVerbose)
      	std::cout << "GBP with " << _gbp.nRegions() << " regions; mem " << _gbp.memory() << "M\n";
      _gbp.setProperties("Schedule=Fixed");
      _gbp.init();
      _gbp.setStopIter(gbpIter); _gbp.setStopObj(gbpObj); _gbp.setStopMsg(-1.0);
      _gbp.setVerbose(doVerbose);
      if (isExact)
      	_gbp.setDamping(-1.0);  // no damping if we think it's exact

      // Get region indices for single-variable beliefs
      mex::vector<mex::gbp::findex> regions(nvar);
      for (size_t v = 0; v<nvar; ++v) {
      	if (!evVar.contains(Var(v, 0))) 
      	  regions[v] = _gbp.regionWith(Var(v, 0));
      }

      double gbpLeft, gbpStop = timeSystem() + gbpTime;
      while ((gbpLeft = gbpStop - timeSystem()) > 0) {
    	_gbp.setStopTime(std::min(dt, gbpLeft));
    	_gbp.run();
    	// save current answer
	switch (task) {
	case Task::PR:{
	  logZ = _gbp.logZStable();
	  res = true;}
	  break;
	case Task::MAR:{
	  for (size_t v = 0; v < nvar; ++v)
	    if (!evVar.contains(Var(v, 0)))
	      bel[v] = _gbp.computeRegionBelief(regions[v]).marginal(Var(v, 0));
	  
	  res = true;}
	  break;
	}
      }
      if (doVerbose)
	std::cout << "GBP " << _gbp.logZ() / c_log10 << "\n";
      if (_gbp.dObj() < gbpObj){
	if (doVerbose)
	  std::cout << "Reached objective tolerance\n";
	res = true; break;
      }
      if (_gbp.iter() >= gbpIter && gbpIter > 0) {
	if (doVerbose)
	  std::cout << "Reached iteration limit\n";
	res = true;break;
      }
      if (_gbp.logZ() == -mex::infty()) {
	if( doVerbose)
	  std::cout << "Model deemed inconsistent\n"; 
	res = false;break;
      }
	
      doneGBP = true;
	
      if (isExact && _gbp.dObj() < gbpObj) { //InducedWidth <= ibound)
	if (doVerbose)
	  std::cout << "Answer should be exact\n"; // TODO:: Continue on to c
	return true;
      }
	
      //TODO;: FIX THIS 
    } catch (std::exception& e) {
      // 	if (_gbp.getDamping() > 0.25) {
      // 	  doneGBP = true;
      // 	  std::cout << "Caught exception (memory overreach?).  Damping on => quitting GBP\n";
      // 	  res = false;
      // 	}
      // 	else {
      // 	  doneGBP = false;
      // 	  MemLimit *= .9; ibound--;
      // 	  std::cout << "Caught exception (memory overreach?).  Trying again with ibound " << ibound << " and MemLimit " << MemLimit << "\n";
      // 	  res = false;
      // 	}
	
    } // end catch



    // }// end try
  } // end try
  return res;

}



bool BpInterface::doIterativeConditioning(){

  if (order.size() == 0)
    order = factGraph.order(mex::graphModel::OrderMethod::MinWidth);  // need an order if none yet...
  mex::gbp _gbp(mex::vector<Factor>());  // !!! blank out GBP object; restore memory
  VarSet cond;

  //!!! TODO: if doCond > 1, we should immediately populate with that many states
  bool isExact = true;
  while (!isExact) {
    // add check for best conditioner given cardinality limit !!!  or, only increment if anytime...
    cond += factGraph.bestConditioner(order, cond);
    if (doVerbose) {
      std::cout << "\n";
      std::cout << "Conditioning " << cond << "\n";
    }
    size_t ibound = iboundInit;																// check all iBounds again (in case higher available)

    bool doneCGBP = false;
    while (!doneCGBP) {
      bool useMBE = false;
      if (task == Task::PR && fitsMBE(factGraph, order, &cond)) {
	if (doVerbose)
	  std::cout << "Trying exact via MBE\n";
	useMBE = true; isExact = true;
      }
      try {
	Factor lnZ(cond);
	bool failed = false;
	mex::vector<mex::vector<Factor> > condMarginals;
	mex::vector<mex::gbp::findex> regions(nvar);			// !!! dangerous; pull out & match to prev (clear+resize)
	for (size_t i = 0; i < lnZ.nrStates(); ++i) {
	  std::map<Var, size_t> val;  ind2sub(cond, i, val);
	  mex::vector<Factor> fcond = factGraph.factors();
	  for (size_t f = 0; f<fcond.size(); ++f) {
	    VarSet isect = cond & fcond[f].vars();
	    if (isect.size() > 0) 
	      fcond[f] = fcond[f].condition(isect, sub2ind(isect, val));
	  }

	  try {
	    if (useMBE) {         // if a bucket elim pass was good enough, do that:
	      mex::graphModel gm(fcond);
	      lnZ[i] = solveMBE(gm, order);
	      for (size_t v = 0; v < cond.size(); ++v) 
		std::cout << cond[v] << "=" << val[cond[v]] << " "; 
	      std::cout << lnZ[i] << "\n";
	      continue;
	    }                     // otherwise we need to do GBP-like updates:
	  }
	  catch (std::exception& e) {
	    std::cout << "Caught exception; failure in MBE; trying GBP\n";
	  }

	  mex::lbp fgcond(fcond);
	  fgcond.setProperties("Schedule=Priority,Distance=L1");
	  fgcond.setStopIter(lbpIter); fgcond.setStopMsg(lbpErr); fgcond.setStopObj(lbpObj);
	  fgcond.setStopTime(0.5);
	  fgcond.init();
	  fgcond.run();
	  fgcond.reparameterize();
	  fcond = fgcond.factors();

	  // TODO:DEBUG: it seems that reinitializing gbp is not equivalent to re-constructing a new one.. FIX
	  mex::gbp _gbp(fcond);
	  if (vm.count("ijgp"))
	    _gbp.setMinimal(true); else _gbp.setMinimal(false);  // use "ijgp-like" regions?
	  isExact = gbpPopulateCliques(_gbp, order, ibound, &cond);
	  _gbp.setProperties("Schedule=Fixed");
	  _gbp.setStopIter(gbpIter); _gbp.setStopObj(gbpObj); _gbp.setStopMsg(-1.0);
	  _gbp.setStopTime(gbpTime); _gbp.setVerbose(doVerbose);
	  _gbp.init();

	  if (task == Task::MAR) {		// TODO: change to loop over "inferred" variable list...
	    for (size_t v = 0; v < nvar; ++v) {
	      if (!evVar.contains(Var(v, 0)) && !cond.contains(Var(v, 0))) 
		regions[v] = _gbp.regionWith(Var(v, 0));
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
	  if (task == Task::MAR) {
	    condMarginals.push_back(bel);
	    for (size_t v = 0; v < nvar; ++v)
	      if (!evVar.contains(Var(v, 0)) && !cond.contains(Var(v, 0)))
		condMarginals[i][v] = _gbp.computeRegionBelief(regions[v]).marginal(Var(v, 0));
	    //condMarginals[i][v]=_gbp.belief(Var(v,0));
	  }
	  for (size_t v = 0; v < cond.size(); ++v) 
	    std::cout << cond[v] << "=" << val[cond[v]] << " "; 
	  std::cout << lnZ[i] << "\n";
	}
	if (failed) { 
	  if (doVerbose)
	    std::cout << "Failing out\n";
	  doneCGBP = true; continue;
	} // if we failed out, condition on more vars

	doneCGBP = true;
	double lnZtot = lnZ.logsumexp();
	switch (task) {
	  //Save solutions
	case Task::PR:
	  logZ = lnZtot;
	  break;
	case Task::MAR:
	  Factor probs = (lnZ - lnZtot).exp();
	  for (size_t v = 0; v < nvar; ++v) {
	    if (evVar.contains(Var(v, 0))) {} // evidence variables not updated
	    else if (cond.contains(Var(v, 0)))  { 
	      bel[v] = probs.marginal(Var(v, 0));
	    }
	    else {		// TODO: change to incremental update?
	      bel[v] = condMarginals[0][v] * probs[0];
	      for (size_t i = 1; i < lnZ.nrStates(); ++i) 
		bel[v] += condMarginals[i][v] * probs[i];
	    }
	  }
	  //writeMAR(outfile, bel);
	  break;
	}
	if (doVerbose)
	  std::cout << "Conditioning " << cond << " => " << lnZtot << " (" << lnZtot / c_log10<< ")\n";

      }
      catch (std::exception& e) {
	doneCGBP = false; MemLimit *= .9; ibound--;
	if (doVerbose)
	  std::cout << "Caught exception (memory overreach?).  Trying again with ibound " << ibound << " and MemLimit " << MemLimit << "\n";
	continue;
	// TODO: this is right if we're only doing this part, but not right if we're running incrementally
      }

      // !!! TODO: if doCond > 1, quit (non-incremental)?
      if (isExact) { 
	if (doVerbose)
	  std::cout << "Answer should be exact\n"; return 0;
      }
    }
  }

  // run CCS expansion on larger regions until timeout (?)

}

/*
  Compute an elimination order
  Params: numTries: the number of orders to compute and score
  timeLimit: max amount of time (in seconds)  to spend searching for an order
  
  returns score of the ordering
*/
double  BpInterface::computeVariableOrder(int numTries, double timeLimit){
  // start with a random order in case the model is very dense
  double score = factGraph.order(mex::graphModel::OrderMethod::Random, order, 0, memUseRandom);

  // Otherwise, calculate elimination order(s) ////////////////////////////////////
  double startOrder = timeSystem();
  size_t iOrder = 0;
  // Try to build new orders until time or count limit reached ////////////////////
  while (iOrder < numTries && (timeSystem() - startOrder < timeLimit)) {
    score = factGraph.order(mex::graphModel::OrderMethod::WtMinFill, order, nExtra, score);
    ++iOrder;
  }
    //if (score < memUseRandom)
  // int InducedWidth = factGraph.inducedWidth(order);
  if(doVerbose)
    std::cout << "Best order of " << iOrder << " has induced width " << factGraph.inducedWidth(order) << ", score " << score << "\n";

  return score;
}


// Utility Functions

bool BpInterface::fitsMBE(const graphModel& gm, const mex::VarOrder& order, const VarSet* cond) {
  double mbCutoff = MemLimit / sizeof(double)* 1024 * 1024;     // translate memory into MBE cutoff
  bool isExact = true;
  mex::mbe mb(gm.factors());
  mb.setOrder(order);
  mb.setProperties("ElimOp=SumUpper,sBound=inf,DoMatch=1,DoMplp=0,DoFill=0,DoJG=0,DoHeur=0");
  mb.setIBound(100); double mbMem = mb.simulateMemory(NULL, cond, mbCutoff, &isExact);

  if (mbMem < mbCutoff && isExact)
    return true;

  return false;
}

double BpInterface::solveMBE(const graphModel& gm, const mex::VarOrder& order) {
  mex::mbe mb(gm.factors());
  mb.setOrder(order);
  mex::vector<mex::mbe::vindex> pt(order.size());
  for (size_t v = 0; v < pt.size(); ++v) pt[v] = -1;
  mb.setPseudotree(pt);
  mb.setProperties("ElimOp=SumUpper,sBound=inf,DoMatch=1,DoMplp=0,DoFill=0,DoJG=0,DoHeur=0");
  mb.setIBound(100); //double mbMem = mb.simulateMemory(NULL,NULL,mbCutoff,&isExact);
  if (doVerbose)
    std::cout << "Attempting exact solve\n";
  //std::cout<<"Attempting exact solve; mbMem="<<mbMem<<" vs "<<mbCutoff<<" ("<<MemLimit<<")\n";
  mb.init();
  return mb.logZ();
}

/**
   Try to compute PR exactly
   
   returns true if successful
 */
bool BpInterface::tryExactPR(const mex::graphModel& gm, const mex::VarOrder& order) {
  try {
    double mbCutoff = MemLimit / sizeof(double)* 1024 * 1024;     // translate memory into MBE cutoff
    isExact = true;
    mex::mbe mb(gm.factors());
    

    ////////////////DELETE ME////////////////////
    printf("MBE Factors\n");
    printFactors(mb._gmo.factors());

    /////////////////////////////////////

    mb.setOrder(order);
    mb.setProperties("ElimOp=SumUpper,sBound=inf,DoMatch=1,DoMplp=0,DoFill=0,DoJG=0,DoHeur=0");
    mb.setIBound(100); double mbMem = mb.simulateMemory(NULL, NULL, mbCutoff, &isExact);
    if (mbMem < mbCutoff && isExact) {
      if (doVerbose)
	std::cout << "Attempting exact solve; mbMem=" << mbMem << " vs " << mbCutoff << " (" << MemLimit << ")\n";
      mb.init();
      logZ = mb.logZ() / c_log10;
      if(doVerbose)
	std::cout << "Exact solution by MBE: " << mb.logZ() / c_log10 << "\n";
      return true;
    }
  }
  catch (std::exception& e) {
    // Failed (probably for memory reasons) => try GBP
    if (doVerbose)
      std::cout << "Failed (due to memory problem?)  Trying GBP\n";
  }
  return false;
}


bool BpInterface::gbpPopulateCliques(mex::gbp& _gbp, const mex::VarOrder& order, size_t& ibound, mex::VarSet* cond) {
  bool isExact = false;
  double mbCutoff = MemLimit / sizeof(double)* 1024 * 1024;     // translate memory into MBE cutoff
  _gbp.clearRegions();
  mex::mbe mb(_gbp.gmOrig().factors());
  mb.setOrder(order);
  mex::vector<mex::mbe::vindex> pt(order.size()); for (size_t v = 0; v < pt.size(); ++v) pt[v] = -1;
  mb.setPseudotree(pt);
  mb.setProperties("ElimOp=SumUpper,sBound=inf,DoMatch=1,DoMplp=0,DoFill=0,DoJG=1,DoHeur=0");
  mb.setIBound(ibound);
  mex::vector<VarSet> cliques;
  double mem = std::numeric_limits<double>::infinity();
  mbCutoff *= 2;  // leave a little slack for mismatch in criteria
  double mbMem = mb.simulateMemory(&cliques, cond, mbCutoff, &isExact);
  if (mbMem < mbCutoff) {
    _gbp.addRegions(cliques); mem = _gbp.memory();
  }
  while (ibound > 0 && (mbMem >= mbCutoff || mem > MemLimit)) {
    if (doVerbose)
      std::cout << "MBE iBound " << ibound << " = " << mem << "M\n";
    mb.setIBound(--ibound); cliques.clear(); mbMem = mb.simulateMemory(&cliques, cond, mbCutoff, &isExact);
    if (mbMem < mbCutoff) {
      _gbp.clearRegions(); _gbp.addRegions(cliques); mem = _gbp.memory();
    }
  }
  //ofstream ofs("cliques.mbe.txt");
  //for (size_t c=0;c<cliques.size();++c) ofs<<cliques[c]<<"\n"; std::cout<<"\n";  // output for DEBUG !!!
  //ofs.close();
  if(doVerbose)
    std::cout << "MBE iBound " << ibound << " = " << mem << "M\n";
  return isExact;
}
