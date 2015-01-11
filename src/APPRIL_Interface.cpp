// VE.cpp : Defines the exported functions for the DLL application.
//

#if defined WINDOWS || _WINDOWS
#include <process.h>    /* _beginthread, _endthread */
#include "stdafx.h"
#endif // WINDOWS

//#define USING_BUCKET_ELIMINATION
//#define USING_DAOOPT

#include "Utils\MiscUtils.hxx"
#include "Utils\Mutex.h"
#include "Problem\Problem.hxx"
#include "CVO\VariableOrderComputation.hxx"

#if defined USING_BUCKET_ELIMINATION
#include <BE\Bucket.hxx>
#include <BE\BEworkspace.hxx>
#endif // USING_BUCKET_ELIMINATION

#if defined USING_DAOOPT
#include "DaooptInterface.h"
#endif // USING_DAOOPT

#include <APPRIL_Interface.hxx>

// this class is used to automatically initialize/terminate the APPRIL interface, e.g. when the module is loaded/unloaded.
class ARRPIL_INTERFACE_AUTO_CONSTR_DESTR
{
public :
	ARRPIL_INTERFACE_AUTO_CONSTR_DESTR(void)
	{
		APPRILinterface::Initialize() ;
	}
	~ARRPIL_INTERFACE_AUTO_CONSTR_DESTR(void)
	{
		APPRILinterface::Terminate() ;
	}
} ;

a single global ARRPIL_INTERFACE_AUTO_CONSTR_DESTR obj; constructor of this obj will initialize the thread, destructor will terminate the thread.
//APPRILinterface::ARRPIL_INTERFACE_AUTO_CONSTR_DESTR APPRIL_GLOBAL_SINGLE_HELPER ;

static long APPRIL_interface_thread_stop_requested = 0 ;
static FILE *APPRIL_fpLOG = NULL ;

#define APPRIL_INTERFACE_CMD_NONE								0
#define APPRIL_INTERFACE_CMD_EXIT								1
#define APPRIL_INTERFACE_CMD_DEFINE_PROBLEM						2
#define APPRIL_INTERFACE_CMD_START_QUERY_COMPUTATION			3
#define APPRIL_INTERFACE_CMD_STOP								4

static const char *MapStateId2String(int cmd_id)
{
	switch (cmd_id) {
		case APPRIL_INTERFACE_STATE_NONE : return "None" ;
		case APPRIL_INTERFACE_STATE_ANALYZING_PROBLEM : return "AnalyzingProblem" ;
		case APPRIL_INTERFACE_STATE_PROBLEM_ANALYZED : return "ProblemAnalyzed" ;
		case APPRIL_INTERFACE_STATE_COMPUTING_QUERY : return "ComputingQuery" ;
		case APPRIL_INTERFACE_STATE_QUERY_COMPUTED : return "QueryComputed" ;
		case APPRIL_INTERFACE_STATE_EXIT_COMPLETED : return "Exiting" ;
		}
	return "" ;
}

static const char *MapCommandId2String(int cmd_id)
{
	switch (cmd_id) {
		case APPRIL_INTERFACE_CMD_NONE : return "None" ;
		case APPRIL_INTERFACE_CMD_EXIT : return "Exit" ;
		case APPRIL_INTERFACE_CMD_DEFINE_PROBLEM : return "DefineProblem" ;
		case APPRIL_INTERFACE_CMD_START_QUERY_COMPUTATION : return "StartQueryComputation" ;
		case APPRIL_INTERFACE_CMD_STOP : return "Stop" ;
		}
	return "" ;
}


const char *APPRILinterface::MapErrorCode2String(int ErrorCode)
{
	switch (ErrorCode) {
		case APPRIL_INTERFACE_ERROR_NONE : return "None" ;
		case APPRIL_INTERFACE_ERROR_UNKNOWN : return "Unknown" ;
		case APPRIL_INTERFACE_ERROR_SYSTEM : return "System" ;
		case APPRIL_INTERFACE_ERROR_NOT_IMPLEMENTED : return "NotImplemented" ;
		case APPRIL_INTERFACE_ERROR_MEMORY_ALLOCATION : return "MemoryAllocation" ;
		case APPRIL_INTERFACE_ERROR_TIMEOUT : return "Timeout" ;
		case APPRIL_INTERFACE_ERROR_FAILED_TO_CREATE_OBJ : return "FailedToCreateObj" ;
		case APPRIL_INTERFACE_ERROR_FAILED_TO_CREATE_CMD : return "FailedToCreateCmd" ;
		case APPRIL_INTERFACE_ERROR_FAILED_TO_POST_CMD : return "FailedToPostCmd" ;
		case APPRIL_INTERFACE_ERROR_PROBLEM_TYPE_NOT_DEFINED : return "ProblemTypeNotDefined" ;
		case APPRIL_INTERFACE_ERROR_PROBLEM_DATA_EMPTY : return "ProblemDataEmpty" ;
		case APPRIL_INTERFACE_ERROR_UNKNOWN_PROBLEM_TYPE : return "UnknownProblemType" ;
		case APPRIL_INTERFACE_ERROR_PROBLEM_UNDEFINED : return "ProblemUndefined" ;
		case APPRIL_INTERFACE_ERROR_FAILED_TO_LOAD_PROBLEM : return "FailedToLoadProblem" ;
		case APPRIL_INTERFACE_ERROR_FAILED_TO_ANALYSE_PROBLEM : return "FailedToAnalyzeProblem" ;
		case APPRIL_INTERFACE_ERROR_PROBLEM_PREPROCESSING : return "FailedToProprocessProblem" ;
		case APPRIL_INTERFACE_ERROR_FAILED_TO_START_CVO : return "FailedToStartVarOrderComputation" ;
		case APPRIL_INTERFACE_ERROR_QUERY_UNDEFINED : return "QueryUndefined" ;
		case APPRIL_INTERFACE_ERROR_QUERY_UNKNOWN : return "QueryUnknown" ;
		case APPRIL_INTERFACE_ERROR_FAILED_TO_CONSTRUCT_SOLUTION_WORKSPACE : return "FailedToConstructSolutionWorkspace" ;
		case APPRIL_INTERFACE_ERROR_MEMORY_BOUND_EXCEEDED : return "MemoryBoundExceeded" ;
		case APPRIL_INTERFACE_ERROR_FAILED_TO_LAUNCH_THREAD : return "FialedToLaunchThread" ;
		case APPRIL_INTERFACE_ERROR_FAILED_TO_LOAD_EVIDENCE : return "FailedToLoadEvidence" ;
		case APPRIL_INTERFACE_ERROR_FAILED_TO_INCORPORATE_EVIDENCE : return "FailedToIncorporateEvidence" ;
		}
	return "NA" ;
}


static int SerializeParameterList(std::list<std::pair<std::string,std::string>> & AssignmentList, char * & BUF, int & lBUF)
{
	std::string s ;
	for (std::list<std::pair<std::string,std::string>>::iterator i = AssignmentList.begin() ; i != AssignmentList.end() ; i++) {
		std::pair<std::string,std::string> & assignment = *i ;
		std::string & sN = assignment.first ;
		std::string & sV = assignment.second ;
		if (s.length() > 0) s += '\n' ;
		s += sN ; s += '=' ; s += sV ;
		}
	if (lBUF < s.length()+1) {
		if (NULL != BUF) { delete [] BUF ; BUF = NULL ; lBUF = 0 ; }
		if (s.length() > 0) { BUF = new char[s.length()+1] ; if (NULL == BUF) return 1 ; lBUF = s.length()+1 ; }
		}
	if (s.length() > 0) {
		memcpy(BUF, s.c_str(), s.length()) ;
		lBUF = s.length() ;
		BUF[lBUF++] = 0 ;
		}
	return 0 ;
}


// ********************************************************************************************************************************
// blackboard
// ********************************************************************************************************************************

class APPRIL_blackboard
{
public :
	ARE::utils::RecursiveMutex _Mutex ;
	int _CurrentState ;
	INT64 _tCurrentState ;
	// engine error; this is used when an error occurs at runtime
	int _CurrentError ;
	// problem data given as input
	std::string _ProblemDataFormat ;
	char *_ProblemData ;
	int _ProblemDataLen ;
	char *_EvidenceData ;
	int _EvidenceDataLen ;
	std::string _FnCombinationOperator ; // e.g. product, sum
	std::string _VarEliminationOperator ; // e.g. sum, min, max
	// ARP problem representation
	ARE::ARP *_Problem ;
	// general parameters
	int _MaxNumProcessorThreads ;
	INT64 _SolutionMemoryBoundInBytes ; // default is 2^30 (1073741824) = 1GB
	// solution complexity
	double _Estimate_NumSearchSpaceNodesExpanded ; // log10 of the actual number; -1 if unknown or NA
	// percentage of query computation work done
	double _SolutionCompletionPercentage ;
	// query answer
	double _Query_Answer ;
	vector<int> _QueryAnswerAssignment ; // e.g. when computing MPE/MAP, the best assignment
	bool _Query_AnswerIsLogScale ;
	int _Query_SingleVariable ;
	int _Query_SingleVariableDomainSize ;
	double _Query_SingleVariableDistribution[32] ;
	// variable order computation
	ARE::VarElimOrderComp::Order _CVObestOrder ; 
	ARE::VarElimOrderComp::CVOcontext _CVOcontext ;
#if defined USING_BUCKET_ELIMINATION
	// bucket elimination
	BucketElimination::BEworkspace _BEws ;
#endif // defined USING_BUCKET_ELIMINATION
#if defined USING_DAOOPT
	daoopt::DaooptInterface _daoopt ;
#endif // defined USING_DAOOPT
public :
	int SetProblemData(const char *Format, const char *Buffer, int Len, const char *EvidenceBuffer, int EvidenceLen)
	{
		ARE::utils::AutoLock lock(_Mutex) ;
		if (NULL != _ProblemData) 
			{ delete [] _ProblemData ; _ProblemData = NULL ; }
		if (NULL != _EvidenceData) 
			{ delete [] _EvidenceData ; _EvidenceData = NULL ; }
		_ProblemDataLen = _EvidenceDataLen = 0 ;
		_ProblemDataFormat.erase() ;
		if (NULL == Buffer || 0 == Len) 
			return 0 ;
		_ProblemData = new char[Len+1] ;
		if (NULL == _ProblemData) 
			return 1 ;
		memcpy(_ProblemData, Buffer, Len) ;
		_ProblemData[Len] = 0 ;
		_ProblemDataLen = Len ;
		if (NULL != Format) 
			_ProblemDataFormat = Format ;
		if (NULL != EvidenceBuffer && EvidenceLen > 0) {
			_EvidenceData = new char[EvidenceLen+1] ;
			if (NULL == _EvidenceData) 
				{ delete [] _ProblemData ; _ProblemData = NULL ; _ProblemDataLen = 0 ; return 1 ; }
			memcpy(_EvidenceData, EvidenceBuffer, EvidenceLen) ;
			_EvidenceData[EvidenceLen] = 0 ;
			_EvidenceDataLen = EvidenceLen ;
			}
		return 0 ;
	}
	int StartQueryComputation(void)
	{
		ARE::utils::AutoLock lock(_Mutex) ;
		_CurrentError = APPRIL_INTERFACE_ERROR_NONE ;
		_SolutionCompletionPercentage = 0.0 ;
		_Query_Answer = DBL_MAX ;
		_QueryAnswerAssignment.clear() ;
		_Query_AnswerIsLogScale = false ;
		_Query_SingleVariableDomainSize = 0 ;
		if (NULL == _Problem) {
			_CurrentError = APPRIL_INTERFACE_ERROR_PROBLEM_UNDEFINED ;
			return APPRIL_INTERFACE_ERROR_FAILED_TO_CONSTRUCT_SOLUTION_WORKSPACE ;
			}
#if defined USING_BUCKET_ELIMINATION
		if (! _BEws.IsValid()) {
			if (NULL != APPRIL_fpLOG) {
				INT64 tNow = ARE::GetTimeInMilliseconds() ;
				fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_bb : error, cannot start BEWS, workspace not set up ...", tNow) ;
				fflush(APPRIL_fpLOG) ;
				}
			_CurrentError = APPRIL_INTERFACE_ERROR_FAILED_TO_CONSTRUCT_SOLUTION_WORKSPACE ;
			return APPRIL_INTERFACE_ERROR_FAILED_TO_CONSTRUCT_SOLUTION_WORKSPACE ;
			}
		if (_SolutionMemoryBoundInBytes >= 0 && _SolutionMemoryBoundInBytes < _BEws.MaxTotalFunctionSpace()) {
			_CurrentError = APPRIL_INTERFACE_ERROR_MEMORY_BOUND_EXCEEDED ;
			return APPRIL_INTERFACE_ERROR_MEMORY_BOUND_EXCEEDED ;
			}
		_BEws.SetFnCombinationType(_Problem->FnCombinationType()) ;
		_BEws.SetVarEliminationType(_Problem->VarEliminationType()) ;
		_BEws.logFile() = APPRIL_fpLOG ;
		_BEws.CreateThread() ;
		if (0 == _BEws._ThreadHandle) {
			if (NULL != APPRIL_fpLOG) {
				INT64 tNow = ARE::GetTimeInMilliseconds() ;
				fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_bb : error, cannot start BEWS, no thread ...", tNow) ;
				fflush(APPRIL_fpLOG) ;
				}
			_CurrentError = APPRIL_INTERFACE_ERROR_FAILED_TO_LAUNCH_THREAD ;
			return APPRIL_INTERFACE_ERROR_FAILED_TO_LAUNCH_THREAD ;
			}
#endif // USING_BUCKET_ELIMINATION
#if defined USING_DAOOPT
		// Nothing to do; when running regular "CheckQueryComputationStatus", a batch of search nodes are expanded.
#endif // USING_DAOOPT
		return 0 ;
	}
	int CheckQueryComputationStatus(bool & IsOngoing)
	{
		IsOngoing = false ;
		ARE::utils::AutoLock lock(_Mutex) ;
		if (NULL == _Problem) 
			return 1 ;
		// check if query computation is complete
#if defined USING_BUCKET_ELIMINATION
		if (0 != _BEws._ThreadHandle) {
			IsOngoing = true ;
			}
		_SolutionCompletionPercentage = _BEws.GetSolutionCompletionPercentage() ;
		_Query_Answer = _BEws.CompleteEliminationResult() ;
		_QueryAnswerAssignment.clear() ;
		_Query_AnswerIsLogScale = false ;
		_Query_SingleVariableDomainSize = _BEws.MarginalSingleVariableDistributionK() ;
		for (int i = 0 ; i < _Query_SingleVariableDomainSize ; i++) 
			_Query_SingleVariableDistribution[i] = _BEws.MarginalSingleVariableDistribution(i) ;
#endif // USING_BUCKET_ELIMINATION
#if defined USING_DAOOPT
		if (APPRIL_INTERFACE_STATE_COMPUTING_QUERY == _CurrentState) {
			bool daoopt_res = _daoopt.solve(1000) ;
			if (daoopt_res) {
				// problem solved
				}
			else {
				IsOngoing = true ;
				}
			}
#endif // USING_DAOOPT
		return 0 ;
	}
	int StopQueryComputation(void)
	{
		ARE::utils::AutoLock lock(_Mutex) ;
		if (NULL == _Problem) 
			return 0 ;
#if defined USING_BUCKET_ELIMINATION
		_BEws.StopThread() ;
		if (NULL != APPRIL_fpLOG) {
			INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
			fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_bb : stop query computation; done ...", tNowLog) ;
			fflush(APPRIL_fpLOG) ;
			}
#endif // USING_BUCKET_ELIMINATION
#if defined USING_DAOOPT
		// Nothing to do
#endif // USING_DAOOPT
		return 0 ;
	}
	int StartProblemAnalysis(void)
	{
		ARE::utils::AutoLock lock(_Mutex) ;
		if (NULL == _Problem) 
			return 1 ;
		if (_Problem->HasVarOrdering()) {
			// no need to launch variable order computation; one is provided as input.
			SetState(APPRIL_INTERFACE_STATE_PROBLEM_ANALYZED) ;
			}
		else {
			_CVOcontext.CreateCVOthread() ;
			if (0 == _CVOcontext._ThreadHandle) {
				if (NULL != APPRIL_fpLOG) {
					INT64 tNow = ARE::GetTimeInMilliseconds() ;
					fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_bb : error, cannot start CVO, no thread ...", tNow) ;
					fflush(APPRIL_fpLOG) ;
					}
				return 1 ;
				}
			SetState(APPRIL_INTERFACE_STATE_ANALYZING_PROBLEM) ;
			}
		return 0 ;
	}
	int CheckProblemAnalysisStatus(bool & IsOngoing)
	{
		IsOngoing = false ;
		ARE::utils::AutoLock lock(_Mutex) ;
		if (NULL == _Problem) 
			return 1 ;
		// check if variable order computation is complete
		if (0 != _CVOcontext._ThreadHandle) {
			IsOngoing = true ;
			}
		return 0 ;
	}
	int StopProblemAnalysis(void)
//	{
//		ARE::utils::AutoLock lock(_Mutex) ;
//		if (NULL == _Problem) 
//			return 0 ;
//		_CVOcontext.StopCVOthread() ;
//		if (NULL != APPRIL_fpLOG) {
//			INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
//			fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_bb : stop variable order computation; done ...", tNowLog) ;
//			fflush(APPRIL_fpLOG) ;
			}
		return 0 ;
	}
	inline int GetState(int & ErrorCode)
	{
		ErrorCode = 0 ;
		ARE::utils::AutoLock lock(_Mutex) ;
		ErrorCode = _CurrentError ;
		return _CurrentState ;
	}
	inline void SetState(int state)
	{
		ARE::utils::AutoLock lock(_Mutex) ;
		_CurrentState = state ;
		_tCurrentState = ARE::GetTimeInMilliseconds() ;
		if (APPRIL_INTERFACE_STATE_PROBLEM_ANALYZED == _CurrentState) {
			// check if problem already has a variable ordering
			if (_Problem->HasVarOrdering()) {
				if (NULL != APPRIL_fpLOG) {
					fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_bb : setting state to \"order-analyzed\", order already has variable order ...", _tCurrentState) ;
					}
				}
			// else copy best order into problem
			else if (_CVObestOrder._Width >= 0 && _CVObestOrder._Width < INT_MAX) {
				if (NULL != APPRIL_fpLOG) {
					fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_bb : setting state to \"order-analyzed\", order w*=%d", _tCurrentState, _CVObestOrder._Width) ;
					}
				int res_copyorder = _Problem->SetVarElimOrdering(_CVObestOrder._VarListInElimOrder, _CVObestOrder._Width) ;
				if (0 != res_copyorder) {
					if (NULL != APPRIL_fpLOG) {
						fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_bb : setting state to \"order-analyzed\", error - failed to copy order to problem, res=%d ...", _tCurrentState, res_copyorder) ;
						goto done ;
						}
					}
				}
			else {
				if (NULL != APPRIL_fpLOG) {
					INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
					fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_bb : error - setting state to \"order-analyzed\", but no var order ..", _tCurrentState) ;
					}
				_CurrentError = APPRIL_INTERFACE_ERROR_FAILED_TO_ANALYSE_PROBLEM ;
				goto done ;
				}


			if (NULL != APPRIL_fpLOG) {
				fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_bb : setting state to \"order-analyzed\", variable order in BT order is :\n    ", _tCurrentState) ;
				for (int i = 0 ; i < _Problem->N() ; i++) 
					fprintf(APPRIL_fpLOG, " %d", (_Problem->VarOrdering_VarList())[i]) ;
				}
#if defined USING_BUCKET_ELIMINATION
			// update BEworkspace; this should generate bucket tree and computation schedule.
			int res_BEwsInit = _BEws.Initialize(*_Problem, NULL, 1) ;
			if (0 != res_BEwsInit) {
				if (NULL != APPRIL_fpLOG) {
					fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_bb : setting state to \"order-analyzed\", error - failed to init BE workspace, res=%d ...", _tCurrentState, res_BEwsInit) ;
					goto done ;
					}
				}
			else {
				if (NULL != APPRIL_fpLOG) {
					fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_bb : setting state to \"order-analyzed\" ok, statistics ...", _tCurrentState) ;
					fprintf(APPRIL_fpLOG, "\n     MaxNumVarsInBucket=%I64d", (INT64) _BEws.MaxNumVarsInBucket()) ;
					fprintf(APPRIL_fpLOG, "\n     TotalOriginalFunctionSpace=%I64d bytes", (INT64) _BEws.TotalOriginalFunctionSpace()) ;
					fprintf(APPRIL_fpLOG, "\n     TotalNewFunctionSpace=%I64d bytes", (INT64) _BEws.TotalNewFunctionSpace()) ;
					fprintf(APPRIL_fpLOG, "\n     TotalNewFunctionComputationComplexity=%I64d operations", (INT64) _BEws.TotalNewFunctionComputationComplexity()) ;
					fprintf(APPRIL_fpLOG, "\n     MaxSimultaneousNewFunctionSpace=%I64d bytes", (INT64) _BEws.MaxSimultaneousNewFunctionSpace()) ;
					fprintf(APPRIL_fpLOG, "\n     MaxTotalFunctionSpace=%I64d bytes", (INT64) _BEws.MaxTotalFunctionSpace()) ;
					}
				}
#endif // USING_BUCKET_ELIMINATION
#if defined USING_DAOOPT
			// Set various options here, as documented in daoopt\ProgramOptions.h
			daoopt::ProgramOptions options ;
			options.ibound = 9;
			options.order_iterations = 25;
			options.seed = 2323;
			options.lds = 1;
			options.problemSpec = _ProblemData ;
			options.problemSpec_len = _ProblemDataLen ;
			options.evidSpec = _EvidenceData ;
			options.evidSpec_len = _EvidenceDataLen ;
			// fetch var elim order
			static std::vector<int> varElimOrder ;
			varElimOrder.clear() ;
			int induced_width = -1 ;
			_Problem->GetVarElimOrdering(varElimOrder, induced_width) ;
			options.varOrder = &varElimOrder ;
			// initialize daoopt interface
			bool res_daoopt_init = _daoopt.initialize() ;
			if (! res_daoopt_init) {
				if (NULL != APPRIL_fpLOG) {
					INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
					fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_bb : error - setting state to \"order-analyzed\"; failed to initialize daoopt ..", _tCurrentState) ;
					}
				_CurrentError = APPRIL_INTERFACE_ERROR_FAILED_TO_CONSTRUCT_SOLUTION_WORKSPACE ;
				goto done ;
				}
			bool res_daoopt_preproc = _daoopt.preprocess(options) ;
			if (! res_daoopt_preproc) {
				if (NULL != APPRIL_fpLOG) {
					INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
					fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_bb : error - setting state to \"order-analyzed\"; failed to preprocess daoopt ..", _tCurrentState) ;
					}
				_CurrentError = APPRIL_INTERFACE_ERROR_FAILED_TO_CONSTRUCT_SOLUTION_WORKSPACE ;
				goto done ;
				}
			_Estimate_NumSearchSpaceNodesExpanded = _daoopt.estimate() ;
			if (0.0 == _Estimate_NumSearchSpaceNodesExpanded) {
				// TODO : problem was solved; save solution
				SetState(APPRIL_INTERFACE_STATE_QUERY_COMPUTED) ;
				_Query_Answer = _daoopt.getSolution(NULL) ;
				_Query_AnswerIsLogScale = true ;
				}
#endif // USING_DAOOPT
			}
		else if (APPRIL_INTERFACE_STATE_QUERY_COMPUTED == _CurrentState) {
			if (NULL != APPRIL_fpLOG) {
				INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
				fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_bb : setting state to \"query-computed\" ok, statistics ...", _tCurrentState) ;
#if defined USING_BUCKET_ELIMINATION
				fprintf(APPRIL_fpLOG, "\n     runtime=%I64d milliseconds", (INT64) _BEws.RunTimeInMilliseconds()) ;
				fprintf(APPRIL_fpLOG, "\n     complete_elimination_result=%g", (double) _BEws.CompleteEliminationResult()) ;
				ARE_Function_TableType *dist = _BEws.MarginalSingleVariableDistribution() ;
				int k = _BEws.MarginalSingleVariableDistributionK() ;
				if (k > 0) {
					fprintf(APPRIL_fpLOG, "\n     marginal_first_var_distribution:v=%d;dist={", (int) _BEws.MarginalSingleVariableDistributionVar()) ;
					for (int j = 0 ; j < k ; j++) {
						if (j > 0) 
							fprintf(APPRIL_fpLOG, ";") ;
						fprintf(APPRIL_fpLOG, "%g", (double) dist[j]) ;
						}
					fprintf(APPRIL_fpLOG, "}") ;
					}
#endif // USING_BUCKET_ELIMINATION
#if defined USING_DAOOPT
				// TODO : problem was solved; save solution
				_Query_Answer = _daoopt.getSolution(&_QueryAnswerAssignment) ;
				_Query_AnswerIsLogScale = true ;
				fprintf(APPRIL_fpLOG, "\n     daoopt_result=%g, assignment={", (double) _Query_Answer) ;
				for (int j = 0 ; j < _QueryAnswerAssignment.size() ; j++) 
					{ if (j > 0) fprintf(APPRIL_fpLOG, ",%d", (int) _QueryAnswerAssignment[j]) ; else fprintf(APPRIL_fpLOG, "%d", (int) _QueryAnswerAssignment[j]) ; }
				fprintf(APPRIL_fpLOG, "}") ;
#endif // USING_DAOOPT
				}
			}
		else {
			// this state has no special processing
			}
done :
		fflush(APPRIL_fpLOG) ;
	}
	int SerializeProblemComplexity(std::list<std::pair<std::string,std::string>> & Answer)
	{
		Answer.clear() ;
		ARE::utils::AutoLock lock(_Mutex) ;
		if (APPRIL_INTERFACE_STATE_NONE == _CurrentState) {
			Answer.push_back(std::pair<std::string,std::string>("error", "problem_undefined")) ;
			return 0 ;
			}
		if (_CVObestOrder._Width >= 0 && _CVObestOrder._Width < INT_MAX && NULL != _CVObestOrder._VarListInElimOrder) {
			char s[64] ;
			sprintf(s, "%d", (int) _CVObestOrder._Width) ;
			Answer.push_back(std::pair<std::string,std::string>("problem_graph_induced_width", s)) ;
			sprintf(s, "%g", (double) _CVObestOrder._Complexity) ;
			Answer.push_back(std::pair<std::string,std::string>("VE_space_complexity_log10", s)) ;
			sprintf(s, "%g", (double) _CVObestOrder._TotalNewFunctionStorageAsNumOfElements) ;
			Answer.push_back(std::pair<std::string,std::string>("VE_new_space_complexity_log10", s)) ;
			}
		return 0 ;
	}
	int SerializeSolutionComplexity(std::list<std::pair<std::string,std::string>> & Answer)
	{
		char s[64] ;
		Answer.clear() ;
		ARE::utils::AutoLock lock(_Mutex) ;
		if (APPRIL_INTERFACE_STATE_NONE == _CurrentState) {
			Answer.push_back(std::pair<std::string,std::string>("error", "problem_undefined")) ;
			return 0 ;
			}
#if defined USING_BUCKET_ELIMINATION
		if (_BEws.IsValid()) {
			sprintf(s, "%d", (int) _BEws.MaxNumVarsInBucket()) ;
			Answer.push_back(std::pair<std::string,std::string>("max_cluster_variable_size", s)) ;
			sprintf(s, "%I64d bytes", (INT64) _BEws.MaxTotalFunctionSpace()) ;
			Answer.push_back(std::pair<std::string,std::string>("space_complexity", s)) ;
			sprintf(s, "%I64d operations", (INT64) _BEws.TotalNewFunctionComputationComplexity()) ;
			Answer.push_back(std::pair<std::string,std::string>("time_complexity", s)) ;
			}
#endif // USING_BUCKET_ELIMINATION
		if (_Estimate_NumSearchSpaceNodesExpanded >= 0.0) {
			sprintf(s, "%g", (double) _Estimate_NumSearchSpaceNodesExpanded) ;
			Answer.push_back(std::pair<std::string,std::string>("num_search_space_nodes_to_expand", s)) ;
			}
		return 0 ;
	}
	int SerializeSolution(std::list<std::pair<std::string,std::string>> & Answer)
	{
		char s[64] ;
		Answer.clear() ;
		ARE::utils::AutoLock lock(_Mutex) ;
		if (APPRIL_INTERFACE_STATE_NONE == _CurrentState) {
			Answer.push_back(std::pair<std::string,std::string>("error", "problem_undefined")) ;
			Answer.push_back(std::pair<std::string,std::string>("query_computation_status", "not_started")) ;
			Answer.push_back(std::pair<std::string,std::string>("query_error", "problem_undefined")) ;
			}
		else if (APPRIL_INTERFACE_STATE_ANALYZING_PROBLEM == _CurrentState) {
			Answer.push_back(std::pair<std::string,std::string>("query_computation_status", "not_started")) ;
			Answer.push_back(std::pair<std::string,std::string>("query_error", "problem_is_being_analyzed")) ;
			}
		else if (APPRIL_INTERFACE_STATE_PROBLEM_ANALYZED == _CurrentState) {
			Answer.push_back(std::pair<std::string,std::string>("query_computation_status", "not_started")) ;
			Answer.push_back(std::pair<std::string,std::string>("query_error", "computation_not_started")) ;
			}
		else if (APPRIL_INTERFACE_STATE_COMPUTING_QUERY != _CurrentState && APPRIL_INTERFACE_STATE_QUERY_COMPUTED != _CurrentState) {
			Answer.push_back(std::pair<std::string,std::string>("query_computation_status", "error")) ;
			Answer.push_back(std::pair<std::string,std::string>("query_error", "unknown")) ;
			}
		else if (APPRIL_INTERFACE_ERROR_NONE != _CurrentError) {
			Answer.push_back(std::pair<std::string,std::string>("query_computation_status", "error")) ;
			Answer.push_back(std::pair<std::string,std::string>("query_error", APPRILinterface::MapErrorCode2String(_CurrentError))) ;
			}
#if defined USING_BUCKET_ELIMINATION
		else if (! _BEws.IsValid()) {
			Answer.push_back(std::pair<std::string,std::string>("query_computation_status", "error")) ;
			Answer.push_back(std::pair<std::string,std::string>("query_error", APPRILinterface::MapErrorCode2String(APPRIL_INTERFACE_ERROR_FAILED_TO_CONSTRUCT_SOLUTION_WORKSPACE))) ;
			}
		else {
			Answer.push_back(std::pair<std::string,std::string>("query_computation_status", NULL != _BEws._ThreadHandle ? "ongoing" : "completed")) ;
			if (_Query_SingleVariable >= 0 && _Query_SingleVariable == _BEws.MarginalSingleVariableDistributionVar()) {
				sprintf(s, "%d", (int) _Query_SingleVariable) ;
				Answer.push_back(std::pair<std::string,std::string>("query_variable", s)) ;
				std::string stemp("{") ;
				for (int i = 0 ; i < _Query_SingleVariableDomainSize ; i++) {
					sprintf(s, "%g", (double) _Query_SingleVariableDistribution[i]) ;
					if (i > 0) 
						stemp += ';' ;
					stemp += s ;
					}
				stemp += "}" ;
				Answer.push_back(std::pair<std::string,std::string>("query_answer", stemp)) ;
				}
			if (_Query_Answer < DBL_MAX) {
				sprintf(s, "%g", (double) _Query_Answer) ;
				if (_Query_AnswerIsLogScale) 
					Answer.push_back(std::pair<std::string,std::string>("partition_function_log10", s)) ;
				else 
					Answer.push_back(std::pair<std::string,std::string>("partition_function", s)) ;
				}
			if (_QueryAnswerAssignment.size() > 0) {
				std::string sAssignment("<") ;
				for (int j = 0 ; j < _QueryAnswerAssignment.size() ; j++) 
					{ if (j > 0) sprintf(s, ",%d", (int) _QueryAnswerAssignment[j]) ; else sprintf(s, "%d", (int) _QueryAnswerAssignment[j]) ; sAssignment += s ; }
				sAssignment += '>' ;
				Answer.push_back(std::pair<std::string,std::string>("query_answer_assignment", sAssignment)) ;
				}
			if (_SolutionCompletionPercentage >= 0.0) {
				sprintf(s, "%g", (double) _SolutionCompletionPercentage) ;
				Answer.push_back(std::pair<std::string,std::string>("query_computation_completion_percentage", s)) ;
				}
			}
#endif // USING_BUCKET_ELIMINATION
#if defined USING_DAOOPT
		else {
			if (_Query_Answer < DBL_MAX) {
				sprintf(s, "%g", (double) _Query_Answer) ;
				if (_Query_AnswerIsLogScale) 
					Answer.push_back(std::pair<std::string,std::string>("query_answer_log10", s)) ;
				else 
					Answer.push_back(std::pair<std::string,std::string>("query_answer", s)) ;
				}
			if (_QueryAnswerAssignment.size() > 0) {
				std::string sAssignment("<") ;
				for (int j = 0 ; j < _QueryAnswerAssignment.size() ; j++) 
					{ if (j > 0) sprintf(s, ",%d", (int) _QueryAnswerAssignment[j]) ; else sprintf(s, "%d", (int) _QueryAnswerAssignment[j]) ; sAssignment += s ; }
				sAssignment += '>' ;
				Answer.push_back(std::pair<std::string,std::string>("query_answer_assignment", sAssignment)) ;
				}
			if (APPRIL_INTERFACE_STATE_QUERY_COMPUTED == _CurrentState) {
				Answer.push_back(std::pair<std::string,std::string>("query_computation_status", "completed")) ;
				Answer.push_back(std::pair<std::string,std::string>("query_computation_completion_percentage", "100")) ;
				}
			else {
				Answer.push_back(std::pair<std::string,std::string>("query_computation_status", "ongoing")) ;
				Answer.push_back(std::pair<std::string,std::string>("query_computation_completion_percentage", "unknown")) ;
				}
			}
#endif // USING_DAOOPT
		return 0 ;
	}
	int Initialize(ARE::ARP *Problem, const std::string & CombOp, const std::string & ElimOp, int QueryVariable)
	{
		ARE::utils::AutoLock lock(_Mutex) ;
		if (NULL != APPRIL_fpLOG) {
			INT64 tNOW = ARE::GetTimeInMilliseconds() ;
			fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_if; Initialize() called ...", tNOW) ;
			fflush(APPRIL_fpLOG) ;
			}
		if (NULL != _Problem) {
			delete _Problem ;
			_Problem = NULL ;
			}
		_Problem = Problem ;
		_Estimate_NumSearchSpaceNodesExpanded = -1.0 ;
		_FnCombinationOperator = CombOp ;
		_VarEliminationOperator = ElimOp ;
		_CurrentError = APPRIL_INTERFACE_ERROR_NONE ;
		_SolutionCompletionPercentage = 0.0 ;
		_Query_Answer = DBL_MAX ;
		_QueryAnswerAssignment.clear() ;
		_Query_AnswerIsLogScale = false ;
		_Query_SingleVariable = QueryVariable ;
		_Query_SingleVariableDomainSize = 0 ;
		// do preprocessing here
		if (NULL == _Problem) {
			_CurrentError = APPRIL_INTERFACE_ERROR_PROBLEM_UNDEFINED ;
			return APPRIL_INTERFACE_ERROR_PROBLEM_UNDEFINED ;
			}
		if (0 == _FnCombinationOperator.length() || 0 == _VarEliminationOperator.length()) {
			delete _Problem ;
			_Problem = NULL ;
			_CurrentError = APPRIL_INTERFACE_ERROR_QUERY_UNDEFINED ;
			return APPRIL_INTERFACE_ERROR_QUERY_UNDEFINED ;
			}
		int res_SetOperators = _Problem->SetOperators(_FnCombinationOperator.c_str(), _VarEliminationOperator.c_str()) ;
		if (0 != res_SetOperators) {
			delete _Problem ;
			_Problem = NULL ;
			_CurrentError = APPRIL_INTERFACE_ERROR_QUERY_UNKNOWN ;
			return APPRIL_INTERFACE_ERROR_QUERY_UNKNOWN ;
			}
		_Problem->SetQueryVariable(QueryVariable) ;
* 2014-03-26 KK : post-construction analysis should be done right after problem is loaded and before evidence is loaded/incorporated.
		int res_InitialAnalysis = _Problem->PerformPostConstructionAnalysis() ;
		if (0 != res_InitialAnalysis) {
			delete _Problem ;
			_Problem = NULL ;
			_CurrentError = APPRIL_INTERFACE_ERROR_PROBLEM_PREPROCESSING ;
			return APPRIL_INTERFACE_ERROR_PROBLEM_PREPROCESSING ;
			}
*/
		int res_ElimSingletonDomVars = _Problem->EliminateSingletonDomainVariables() ;
		if (0 != res_ElimSingletonDomVars) {
			delete _Problem ;
			_Problem = NULL ;
			_CurrentError = APPRIL_INTERFACE_ERROR_PROBLEM_PREPROCESSING ;
			return APPRIL_INTERFACE_ERROR_PROBLEM_PREPROCESSING ;
			}
		// initialize CVO stuff
		_CVOcontext._fpLOG = APPRIL_fpLOG ;
		_CVOcontext._Problem = _Problem ;
		_CVOcontext._FindPracticalVariableOrder = true ;
		_CVOcontext._AlgCode = ARE::VarElimOrderComp::MinFill ;
		_CVOcontext._ObjCode = ARE::VarElimOrderComp::Width ;
		_CVOcontext._nRunsToDoMin = 100 ;
		_CVOcontext._nRunsToDoMax = 1000 ;
		_CVOcontext._TimeLimitInMilliSeconds = 60000 ;
		_CVOcontext._nRandomPick = 8 ;
		_CVOcontext._eRandomPick = 0.5 ;
		_CVOcontext._EarlyTerminationOfBasic_W = true ;
		_CVOcontext._EarlyTerminationOfBasic_C = false ;
		_CVOcontext._nThreads = _MaxNumProcessorThreads > 1 ? _MaxNumProcessorThreads-1 : 1 ;
		leave BE workspace as is (should be in reset state); when var order is computed, BEws should be initialized.
		int res_startcvo = StartProblemAnalysis() ;
		if (0 != res_startcvo) {
			delete _Problem ;
			_Problem = NULL ;
			_CurrentError = APPRIL_INTERFACE_ERROR_FAILED_TO_START_CVO ;
			return APPRIL_INTERFACE_ERROR_FAILED_TO_START_CVO ;
			}
		return 0 ;
	}
	int Destroy(void)
	{
		ARE::utils::AutoLock lock(_Mutex) ;
		reset CVOcontext/CVObestorder
		_CVOcontext.Destroy() ;
		_CVObestOrder.Initialize(0) ;
#if defined USING_BUCKET_ELIMINATION
		reset BEworkspace
		_BEws.Destroy() ;
#endif // USING_BUCKET_ELIMINATION
		reset locally
		if (NULL != _Problem) {
			delete _Problem ;
			_Problem = NULL ;
			}
		if (NULL != _ProblemData) 
			{ delete [] _ProblemData ; _ProblemData = NULL ; }
		if (NULL != _EvidenceData) 
			{ delete [] _EvidenceData ; _EvidenceData = NULL ; }
		_ProblemDataLen = _EvidenceDataLen = 0 ;
		_ProblemDataFormat.erase() ;
		return 0 ;
	}
public :
	APPRIL_blackboard(void) 
		:
		_CurrentState(APPRIL_INTERFACE_STATE_NONE), 
		_tCurrentState(0), 
		_ProblemData(NULL), 
		_ProblemDataLen(0), 
		_EvidenceData(NULL), 
		_EvidenceDataLen(0), 
		_Problem(NULL), 
		_MaxNumProcessorThreads(1), 
		_SolutionMemoryBoundInBytes(1073741824), 
		_Estimate_NumSearchSpaceNodesExpanded(-1.0) 
#if defined USING_BUCKET_ELIMINATION
		, _BEws(NULL) 
#endif // USING_BUCKET_ELIMINATION
	{
		CVO setup
		_CVOcontext._BestOrder = &_CVObestOrder ;
		general setup
#if defined _WINDOWS || WINDOWS
		{
			SYSTEM_INFO sysinfo ;
			GetSystemInfo(&sysinfo) ;
			_MaxNumProcessorThreads = sysinfo.dwNumberOfProcessors ;
		}
#endif
#ifdef _DEBUG
		_MaxNumProcessorThreads = 1 ;
#endif // _DEBUG
	}
	~APPRIL_blackboard(void) 
	{
		Destroy() ;
	}
} ;

static APPRIL_blackboard APPRIL_bb ;

********************************************************************************************************************************
Command message stuff
********************************************************************************************************************************

static ARE::utils::RecursiveMutex APPRIL_CMD_Mutex ;

class APPRIL_cmd_msg
{
public :
	int _CMD ;
	ARE::ARP *_Problem ;
	std::string _FnCombinationOperator ; // e.g. product, sum
	std::string _VarEliminationOperator ; // e.g. sum, min, max
	int _QueryVariable ;
	APPRIL_cmd_msg *_Next ;
	parameters
	INT64 _SolutionMemoryBoundInBytes ;
public :
	APPRIL_cmd_msg(int CMD = APPRIL_INTERFACE_CMD_NONE) 
		:
		_CMD(CMD), 
		_Problem(NULL), 
		_QueryVariable(-1), 
		_Next(NULL), 
		_SolutionMemoryBoundInBytes(-1)
	{
	}
	~APPRIL_cmd_msg(void) 
	{
		if (NULL != _Problem) 
			delete _Problem ;
	}
} ;

static APPRIL_cmd_msg *CMD_MSG_QUEUE_HEAD = NULL, *CMD_MSG_QUEUE_TAIL = NULL ;

static int PostCmdMsg(APPRIL_cmd_msg & CMD, bool DumpCurrentQueue)
{
	if (0 != APPRIL_interface_thread_stop_requested) {
		return 1 ;
		}

	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_if; PostCmdMsg id=%s ...", tNOW, MapCommandId2String(CMD._CMD)) ;
		fflush(APPRIL_fpLOG) ;
		}

	CMD._Next = NULL ; // safety

	switch (CMD._CMD) { // check cmd code is valid
		case APPRIL_INTERFACE_CMD_EXIT : // don't post this cmd; signal exit using variable.
			APPRIL_interface_thread_stop_requested = 1 ;
			delete &CMD ;
			return 0 ;
		case APPRIL_INTERFACE_CMD_DEFINE_PROBLEM :
		case APPRIL_INTERFACE_CMD_START_QUERY_COMPUTATION :
		case APPRIL_INTERFACE_CMD_STOP :
			break ;
		default :
			return 1 ;
		}

	{
	ARE::utils::AutoLock lock(APPRIL_CMD_Mutex) ;

	if (DumpCurrentQueue) {
		while (NULL != CMD_MSG_QUEUE_HEAD) {
			APPRIL_cmd_msg *cmd = CMD_MSG_QUEUE_HEAD->_Next ;
			CMD_MSG_QUEUE_HEAD->_Next = NULL ;
			delete CMD_MSG_QUEUE_HEAD ;
			CMD_MSG_QUEUE_HEAD = cmd ;
			}
		CMD_MSG_QUEUE_HEAD = CMD_MSG_QUEUE_TAIL = NULL ;
		}

	if (NULL != CMD_MSG_QUEUE_TAIL) {
		CMD_MSG_QUEUE_TAIL->_Next = &CMD ;
		CMD_MSG_QUEUE_TAIL = &CMD ;
		}
	else {
		CMD_MSG_QUEUE_HEAD = CMD_MSG_QUEUE_TAIL = &CMD ;
		}
	}

	return 0 ;
}


static APPRIL_cmd_msg *PopCmdMsg(APPRIL_cmd_msg *Pop_Only_iff_CMD_Is_This)
{
	APPRIL_cmd_msg *cmd = NULL ;
	ARE::utils::AutoLock lock(APPRIL_CMD_Mutex) ;
	if (NULL != CMD_MSG_QUEUE_HEAD) {
		if (NULL != Pop_Only_iff_CMD_Is_This) {
			if (Pop_Only_iff_CMD_Is_This != CMD_MSG_QUEUE_HEAD) 
				return NULL ;
			}
		cmd = CMD_MSG_QUEUE_HEAD ;
		CMD_MSG_QUEUE_HEAD = CMD_MSG_QUEUE_HEAD->_Next ;
		if (NULL == CMD_MSG_QUEUE_HEAD) 
			CMD_MSG_QUEUE_TAIL = NULL ;
		cmd->_Next = NULL ;
		}
	return cmd ;
}


static int EmptyCmdMsgQueue(void)
{
	ARE::utils::AutoLock lock(APPRIL_CMD_Mutex) ;
	while (NULL != CMD_MSG_QUEUE_HEAD) {
		APPRIL_cmd_msg *cmd = CMD_MSG_QUEUE_HEAD->_Next ;
		CMD_MSG_QUEUE_HEAD->_Next = NULL ;
		delete CMD_MSG_QUEUE_HEAD ;
		CMD_MSG_QUEUE_HEAD = cmd ;
		}
	CMD_MSG_QUEUE_HEAD = CMD_MSG_QUEUE_TAIL = NULL ;
	return 0 ;
}


********************************************************************************************************************************
APPRIL interface thread
********************************************************************************************************************************

#if defined WINDOWS || _WINDOWS
static uintptr_t APPRIL_interface_ThreadHandle = 0 ;
#elif defined (LINUX)
static pthread_t APPRIL_interface_ThreadHandle = 0 ;
#endif 

#if defined WINDOWS || _WINDOWS
typedef unsigned int (*pAPPRIL_Interface_ThreadFn)(void *X) ;
static unsigned int __stdcall APPRIL_Interface_ThreadFn(void *X) 
#elif defined (LINUX)
typedef void *(*pAPPRIL_Interface_ThreadFn)(void *X) ;
static void *APPRIL_Interface_ThreadFn(void *X)
#endif 
{
	if (NULL != APPRIL_fpLOG) {
		INT64 tNow = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : start ...", tNow) ;
		fflush(APPRIL_fpLOG) ;
		}

	/*
		cmd description :

		"exit" - stop all; destroy data; exit.
		"stop" - stop current activity; keep data alive.
		"define problem" - a new problem is being defined; stop current activity, destroy data, accept new problem, start pre-processing.
		"start query" - compute the query for the current problem. this may fail is query cannot be computed.

		states :

		NONE = no current problem; blackboard is blank; thread is idle.
		ANALYZING_PROBLEM = computing variable elimination ordering. thread goes to this state when a new problem is defined. order data is posted on the blackboard.
		PROBLEM_ANALYZED = variable elimination ordering is computed. we are ready to start query computation. order data is on blackboard.
		COMPUTING_QUERY = query is being computed. any results (intermediate) are on the blackboard.
		QUERY_COMPUTED = query is computed. results are on the blackboard.
		EXIT_COMPLETED = thread has exited (or is about to); this state should only be seen at the very end.

		target state upon receipt of cmd :

		"exit" 
			current state = any -> NONE
		"define"
			current state = any -> NONE
		"stop"
			current state = ANALYZING_PROBLEM -> PROBLEM_ANALYZED
			current state = COMPUTING_QUERY -> QUERY_COMPUTED
			otherwise -> maintain current state
		"start" 
			current state = ANALYZING_PROBLEM -> PROBLEM_ANALYZED
			otherwise -> maintain current state

		state transition : (current state X target state); we will ignore cases when we should do nothing.

		EXIT_COMPLETED
			target state ignored; thread will exit.
		NONE 
			EXIT_COMPLETED - switch state; exit thread.
		ANALYZING_PROBLEM 
			NONE - stop computation & reset; when done, switch state.
			PROBLEM_ANALYZED - stop computation, keep current data (blackboard); when done, switch state.
		PROBLEM_ANALYZED 
			NONE - reset; when done, switch state.
			COMPUTING_QUERY - start query; switch state.
		COMPUTING_QUERY 
			NONE - stop computation & reset; when done, switch state.
			QUERY_COMPUTED - stop computation; when done switch state.
		QUERY_COMPUTED 
			NONE - reset; when done, switch state.

		basic logic : 

			when there is cmd, take a peek. based on the cmd id, generate target state; the idea here being that when the current_state=target_state, 
			the cmd can be execured. if no next command, target state is NA.
			compare current state and target state; if different, take appropriate action to get to target state.
			if current_state=target_state, pop cmd and execute.
	*/

	int target_state = APPRIL_INTERFACE_STATE_NONE, next_cmd_id = APPRIL_INTERFACE_CMD_NONE ;
	APPRIL_cmd_msg *cmd2execute = NULL, *next_cmd = NULL ;
	INT64 tNow = 0 ;
	bool stop_exit_logged = false ;

	int & current_state = APPRIL_bb._CurrentState ; current_state = APPRIL_INTERFACE_STATE_NONE ;
	INT64 & tCurrentState = APPRIL_bb._tCurrentState ; tCurrentState = 0 ;

	try {
	while (true) {
		target_state = -1 ;
		next_cmd_id = -1 ;
		cmd2execute = NULL ;
		next_cmd = NULL ;
		if (0 != APPRIL_interface_thread_stop_requested) {
			if (NULL != APPRIL_fpLOG && ! stop_exit_logged) {
				stop_exit_logged = true ;
				INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
				fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : stop/exit requested ...", tNowLog) ;
				fflush(APPRIL_fpLOG) ;
				}
			target_state = APPRIL_INTERFACE_STATE_NONE ;
			if (APPRIL_INTERFACE_STATE_NONE == current_state) 
				goto done ;
			}
		else {
			ARE::utils::AutoLock lock(APPRIL_CMD_Mutex) ;
			peek at the next cmd
			if (NULL != CMD_MSG_QUEUE_HEAD) {
				next_cmd = CMD_MSG_QUEUE_HEAD ;
				next_cmd_id = CMD_MSG_QUEUE_HEAD->_CMD ;
				switch (CMD_MSG_QUEUE_HEAD->_CMD) {
					case APPRIL_INTERFACE_CMD_DEFINE_PROBLEM :
						target_state = APPRIL_INTERFACE_STATE_NONE ;
						break ;
					case APPRIL_INTERFACE_CMD_STOP :
						if (APPRIL_INTERFACE_STATE_ANALYZING_PROBLEM == current_state) 
							we can start executing the cmd right away
							target_state = current_state ;
						else if (APPRIL_INTERFACE_STATE_COMPUTING_QUERY == current_state) 
							we can start executing the cmd right away
							target_state = current_state ;
						else cannot execute the cmd : either query already computed or incompatible state (e.g. problem not defined)
						break ;
					case APPRIL_INTERFACE_CMD_START_QUERY_COMPUTATION :
						if (APPRIL_INTERFACE_STATE_ANALYZING_PROBLEM == current_state) 
							need to stop order computation first, only then can execute the cmd
							target_state = APPRIL_INTERFACE_STATE_PROBLEM_ANALYZED ;
						else if (APPRIL_INTERFACE_STATE_PROBLEM_ANALYZED == current_state) 
							we can start executing the cmd right away
							target_state = current_state ;
						else cannot execute the cmd : either query already computed or incompatible state (e.g. problem not defined)
						break ;
					case APPRIL_INTERFACE_CMD_EXIT :
						target_state = APPRIL_INTERFACE_STATE_NONE ;
						break ;
					case APPRIL_INTERFACE_CMD_NONE : // this is not supposed to be used
					default : // unknown command, or we forgot to handle some command.
						break ;
					}
				if target_state is UNDEFINED, the cmd is invalid in current state; pop/delete next cmd.
				if (target_state < 0) {
					cmd2execute = CMD_MSG_QUEUE_HEAD ;
					log invalid command
					if (NULL != APPRIL_fpLOG) {
						INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
						fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : invalid command; current_state=%s, next_cmd=%s", tNowLog, MapStateId2String(current_state), MapCommandId2String(cmd2execute->_CMD)) ;
						fflush(APPRIL_fpLOG) ;
						}
					CMD_MSG_QUEUE_HEAD = CMD_MSG_QUEUE_HEAD->_Next ;
					if (NULL == CMD_MSG_QUEUE_HEAD) 
						CMD_MSG_QUEUE_TAIL = NULL ;
					cmd2execute->_Next = NULL ;
					delete cmd2execute ;
					continue ;
					}
				if target_state == current_state, the next cmd can be executed; pop it here.
				if (current_state == target_state) {
					log command popped/ready for execution
					cmd2execute = CMD_MSG_QUEUE_HEAD ;
					if (NULL != APPRIL_fpLOG) {
						INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
						fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : ready to execute command; current_state=%s, next_cmd=%s", tNowLog, MapStateId2String(current_state), MapCommandId2String(cmd2execute->_CMD)) ;
						fflush(APPRIL_fpLOG) ;
						}
					CMD_MSG_QUEUE_HEAD = CMD_MSG_QUEUE_HEAD->_Next ;
					if (NULL == CMD_MSG_QUEUE_HEAD) 
						CMD_MSG_QUEUE_TAIL = NULL ;
					cmd2execute->_Next = NULL ;
					}
				}
			}

		if command was popped, execute it. we assume execution is possible (current_state == target_state).
		handle only commands that require particular 3rd party action.
		if (NULL != cmd2execute) {
			switch (cmd2execute->_CMD) {
				case APPRIL_INTERFACE_CMD_DEFINE_PROBLEM :
					{
					if (NULL != APPRIL_fpLOG) {
						INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
						fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : executing \"define\"; current_state=%s", tNowLog, MapStateId2String(current_state)) ;
						fflush(APPRIL_fpLOG) ;
						}
					int res_init = APPRIL_bb.Initialize(cmd2execute->_Problem, cmd2execute->_FnCombinationOperator, cmd2execute->_VarEliminationOperator, cmd2execute->_QueryVariable) ;
					if (0 != res_init) {
						if (NULL != APPRIL_fpLOG) {
							INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
							fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : executing \"define\"; current_state=%s; BB_init failed, res=%s ...", tNowLog, MapStateId2String(current_state), APPRILinterface::MapErrorCode2String(res_init)) ;
							fflush(APPRIL_fpLOG) ;
							}
						APPRIL_bb.SetState(APPRIL_INTERFACE_STATE_NONE) ;
						}
					else {
						if (NULL != APPRIL_fpLOG) {
							INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
							fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : executing \"define\"; current_state=%s; succeeded ...", tNowLog, MapStateId2String(current_state)) ;
							fflush(APPRIL_fpLOG) ;
							}
						}
					cmd2execute->_Problem = NULL ;
					}
					break ;
				case APPRIL_INTERFACE_CMD_STOP :
					if (APPRIL_INTERFACE_STATE_ANALYZING_PROBLEM == current_state) {
						stop variable elimination order computation; keep BB. switch state when done.
						if (NULL != APPRIL_fpLOG) {
							INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
							fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : executing \"stop-cvo\"; current_state=%s", tNowLog, MapStateId2String(current_state)) ;
							fflush(APPRIL_fpLOG) ;
							}
						APPRIL_bb.StopProblemAnalysis() ;
						APPRIL_bb.SetState(APPRIL_INTERFACE_STATE_PROBLEM_ANALYZED) ;
						}
					else if (APPRIL_INTERFACE_STATE_COMPUTING_QUERY == current_state) {
						stop query computation; keep BB. switch state when done.
						if (NULL != APPRIL_fpLOG) {
							INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
							fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : executing \"stop-query\"; current_state=%s", tNowLog, MapStateId2String(current_state)) ;
							fflush(APPRIL_fpLOG) ;
							}
						APPRIL_bb.StopQueryComputation() ;
						APPRIL_bb.SetState(APPRIL_INTERFACE_STATE_QUERY_COMPUTED) ;
						}
					else {
						cannot stop, current state incompatible.
						if (NULL != APPRIL_fpLOG) {
							INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
							fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : executing \"stop\"; current_state=%s invalid, will ignore ...", tNowLog, MapStateId2String(current_state)) ;
							fflush(APPRIL_fpLOG) ;
							}
						}
					break ;
				case APPRIL_INTERFACE_CMD_START_QUERY_COMPUTATION :
					if (APPRIL_INTERFACE_STATE_PROBLEM_ANALYZED != current_state) {
						if (NULL != APPRIL_fpLOG) {
							INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
							fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : executing \"startquery\"; current_state=%s invalid, will ignore ...", tNowLog, MapStateId2String(current_state)) ;
							fflush(APPRIL_fpLOG) ;
							}
						}
					else {
						if (NULL != APPRIL_fpLOG) {
							INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
							fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : executing \"startquery\"; current_state=%s, next_cmd=%s", tNowLog, MapStateId2String(current_state), MapCommandId2String(cmd2execute->_CMD)) ;
							fflush(APPRIL_fpLOG) ;
							}
						if (cmd2execute->_SolutionMemoryBoundInBytes >= 0) 
							APPRIL_bb._SolutionMemoryBoundInBytes = cmd2execute->_SolutionMemoryBoundInBytes ;
						int res = APPRIL_bb.StartQueryComputation() ;
						if (0 != res) {
							if (NULL != APPRIL_fpLOG) {
								INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
								fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : executing \"startquery\"; current_state=%s; failed, res=%d ...", tNowLog, MapStateId2String(current_state), APPRILinterface::MapErrorCode2String(res)) ;
								fflush(APPRIL_fpLOG) ;
								}
							query is computed, but with error code.
							APPRIL_bb.SetState(APPRIL_INTERFACE_STATE_QUERY_COMPUTED) ;
							}
						else {
							APPRIL_bb.SetState(APPRIL_INTERFACE_STATE_COMPUTING_QUERY) ;
							}
						}
					break ;
				default :
					break ;
				}
			delete cmd2execute ;
			goto kill_some_time ;
			}
		otherwise check if action is needed so that preconditions of next cmd are satisfied. we can assume that current_state != target_state.
		else if (target_state >= 0) {
			int have_action = -1 ;
			switch (current_state) {
				case APPRIL_INTERFACE_STATE_ANALYZING_PROBLEM :
					switch (target_state) {
						case APPRIL_INTERFACE_STATE_NONE :
							stop variable elimination order computation; delete BB. switch state when done.
							if (NULL != APPRIL_fpLOG) {
								INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
								fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : action : current_state=%s, target_state=%s", tNowLog, MapStateId2String(current_state), MapStateId2String(target_state)) ;
								fflush(APPRIL_fpLOG) ;
								}
							APPRIL_bb.StopProblemAnalysis() ;
							APPRIL_bb.Destroy() ;
							APPRIL_bb.SetState(APPRIL_INTERFACE_STATE_NONE) ;
							have_action = 1 ;
							break ;
						case APPRIL_INTERFACE_STATE_PROBLEM_ANALYZED :
							stop variable elimination order computation; keep BB. switch state when done.
							if (NULL != APPRIL_fpLOG) {
								INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
								fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : action : current_state=%s, target_state=%s", tNowLog, MapStateId2String(current_state), MapStateId2String(target_state)) ;
								fflush(APPRIL_fpLOG) ;
								}
							APPRIL_bb.StopProblemAnalysis() ;
							APPRIL_bb.SetState(APPRIL_INTERFACE_STATE_PROBLEM_ANALYZED) ;
							have_action = 1 ;
							break ;
						default :
							have_action = 0 ;
							break ;
						}
					break ;
				case APPRIL_INTERFACE_STATE_PROBLEM_ANALYZED :
					switch (target_state) {
						case APPRIL_INTERFACE_STATE_NONE :
							delete BB. switch state when done.
							if (NULL != APPRIL_fpLOG) {
								INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
								fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : action : current_state=%s, target_state=%s", tNowLog, MapStateId2String(current_state), MapStateId2String(target_state)) ;
								fflush(APPRIL_fpLOG) ;
								}
							APPRIL_bb.Destroy() ;
							APPRIL_bb.SetState(APPRIL_INTERFACE_STATE_NONE) ;
							have_action = 1 ;
							break ;
						case APPRIL_INTERFACE_STATE_COMPUTING_QUERY :
							start query computation. switch state when done.
							if (NULL != APPRIL_fpLOG) {
								INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
								fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : action : current_state=%s, target_state=%s", tNowLog, MapStateId2String(current_state), MapStateId2String(target_state)) ;
								fflush(APPRIL_fpLOG) ;
								}
							APPRIL_bb.StartQueryComputation() ;
							APPRIL_bb.SetState(APPRIL_INTERFACE_STATE_COMPUTING_QUERY) ;
							have_action = 1 ;
							break ;
						default :
							have_action = 0 ;
							break ;
						}
					break ;
				case APPRIL_INTERFACE_STATE_COMPUTING_QUERY :
					switch (target_state) {
						case APPRIL_INTERFACE_STATE_NONE :
							stop query computation; delete BB. switch state when done.
							if (NULL != APPRIL_fpLOG) {
								INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
								fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : action : current_state=%s, target_state=%s", tNowLog, MapStateId2String(current_state), MapStateId2String(target_state)) ;
								fflush(APPRIL_fpLOG) ;
								}
							APPRIL_bb.StopQueryComputation() ;
							APPRIL_bb.Destroy() ;
							APPRIL_bb.SetState(APPRIL_INTERFACE_STATE_NONE) ;
							have_action = 1 ;
							break ;
						case APPRIL_INTERFACE_STATE_QUERY_COMPUTED :
							stop query computation; keep BB. switch state when done.
							if (NULL != APPRIL_fpLOG) {
								INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
								fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : action : current_state=%s, target_state=%s", tNowLog, MapStateId2String(current_state), MapStateId2String(target_state)) ;
								fflush(APPRIL_fpLOG) ;
								}
							APPRIL_bb.StopQueryComputation() ;
							APPRIL_bb.SetState(APPRIL_INTERFACE_STATE_QUERY_COMPUTED) ;
							have_action = 1 ;
							break ;
						default :
							have_action = 0 ;
							break ;
						}
					break ;
				case APPRIL_INTERFACE_STATE_QUERY_COMPUTED :
					switch (target_state) {
						case APPRIL_INTERFACE_STATE_NONE :
							delete BB. switch state when done.
							if (NULL != APPRIL_fpLOG) {
								INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
								fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : action : current_state=%s, target_state=%s", tNowLog, MapStateId2String(current_state), MapStateId2String(target_state)) ;
								fflush(APPRIL_fpLOG) ;
								}
							APPRIL_bb.Destroy() ;
							APPRIL_bb.SetState(APPRIL_INTERFACE_STATE_NONE) ;
							have_action = 1 ;
							break ;
						default :
							have_action = 0 ;
							break ;
						}
					break ;
				case APPRIL_INTERFACE_STATE_EXIT_COMPLETED :
				case APPRIL_INTERFACE_STATE_NONE :
					have_action = 0 ;
					break ;
				}
			if (0 == have_action) {
				if (NULL != APPRIL_fpLOG) {
					INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
					fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : no valid action for state transition; current_state=%d, target_state=%d", tNowLog, (int) current_state, (int) target_state) ;
					fflush(APPRIL_fpLOG) ;
					}
				}
			else if (have_action < 0) {
				if (NULL != APPRIL_fpLOG) {
					INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
					fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : undefined action for state transition; current_state=%d, target_state=%d", tNowLog, (int) current_state, (int) target_state) ;
					fflush(APPRIL_fpLOG) ;
					}
				}
			if (have_action <= 0) {
				we have no action to reach the target_state. delete the cmd.
				APPRIL_cmd_msg *cmd = PopCmdMsg(next_cmd) ;
				if (NULL != cmd) 
					delete cmd ;
				}
			else {
				if (NULL != APPRIL_fpLOG) {
					INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
					fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : action executed : current_state=%s", tNowLog, MapStateId2String(current_state)) ;
					fflush(APPRIL_fpLOG) ;
					}
				hopefully target_state == current_state so that we can execute the cmd
				goto run_again ;
				}
			}

		check if current activity has finished and we can change state
		switch (current_state) {
			bool is_ongoing ; is_ongoing = false ;
			case APPRIL_INTERFACE_STATE_ANALYZING_PROBLEM :
				check if computation has finished, update state.
				APPRIL_bb.CheckProblemAnalysisStatus(is_ongoing) ;
				if (! is_ongoing) {
					APPRIL_bb.SetState(APPRIL_INTERFACE_STATE_PROBLEM_ANALYZED) ;
					}
				break ;
			case APPRIL_INTERFACE_STATE_COMPUTING_QUERY :
				check if query computation has finished, update state.
				APPRIL_bb.CheckQueryComputationStatus(is_ongoing) ;
				if (! is_ongoing) {
					APPRIL_bb.SetState(APPRIL_INTERFACE_STATE_QUERY_COMPUTED) ;
					}
				break ;
			case APPRIL_INTERFACE_STATE_PROBLEM_ANALYZED :
			case APPRIL_INTERFACE_STATE_QUERY_COMPUTED :
			case APPRIL_INTERFACE_STATE_EXIT_COMPLETED :
			case APPRIL_INTERFACE_STATE_NONE :
				do nothing.
				break ;
			}

kill_some_time :
		Sleep(50) ;
run_again : ;
		}
		} catch (...) {
			if (NULL != APPRIL_fpLOG) {
				INT64 tNow = ARE::GetTimeInMilliseconds() ;
				fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : exception ...", tNow) ;
				fflush(APPRIL_fpLOG) ;
				}
			}

done :
	if (NULL != APPRIL_fpLOG) {
		INT64 tNow = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL_th : end ...", tNow) ;
		fflush(APPRIL_fpLOG) ;
		}
	current_state = APPRIL_INTERFACE_STATE_EXIT_COMPLETED ;
	APPRIL_interface_ThreadHandle = 0 ;
#if defined WINDOWS || _WINDOWS
	_endthreadex(0) ;
	return 0  ;
#else
	return NULL ;
#endif
}

int APPRILinterface::Initialize(bool logStuff)
{
	if (logStuff) {
		APPRIL_fpLOG = fopen("APPRIL_log.txt", "w") ;
		}
	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(Initialize); version=%s, build date/time %s %s ...", tNOW, APPRIL_INTERFACE_VERSION_STRING, __DATE__, __TIME__) ;
		fflush(APPRIL_fpLOG) ;
		}

	if (0 != APPRIL_interface_ThreadHandle) 
		return 0 ;

	APPRIL_interface_thread_stop_requested = 0 ;
#if defined WINDOWS || _WINDOWS
	APPRIL_interface_ThreadHandle = _beginthreadex(NULL, 0, APPRIL_Interface_ThreadFn, NULL, 0, NULL) ;
#else
	pthread_create(&APPRIL_interface_ThreadHandle, NULL, APPRIL_Interface_ThreadFn, NULL) ; // TODO third argument
#endif 
	if (0 == APPRIL_interface_ThreadHandle) 
		return APPRIL_INTERFACE_ERROR_SYSTEM ;

	return APPRIL_INTERFACE_ERROR_NONE ;
}


int APPRILinterface::Terminate(void)
{
	APPRIL_interface_thread_stop_requested = 1 ;

	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(Terminate) ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	if (NULL == APPRIL_interface_ThreadHandle) 
		return 0 ;

	INT64 tB = ARE::GetTimeInMilliseconds() ;
	bool hard_kill = false ;
	while (NULL != APPRIL_interface_ThreadHandle) {
		INT64 tNow = ARE::GetTimeInMilliseconds() ;
		INT64 dt = tNow - tB ;
#ifdef _DEBUG
		if (dt > 10000000) {
#else
		if (dt > 10000) {
#endif // _DEBUG
			ok, waited long enough; do hard kill.
			hard_kill = true ;
			TerminateThread((HANDLE) APPRIL_interface_ThreadHandle, APPRIL_INTERFACE_ERROR_TIMEOUT) ;
			CloseHandle((HANDLE) APPRIL_interface_ThreadHandle) ;
			if (NULL != APPRIL_fpLOG) {
				INT64 tNow = ARE::GetTimeInMilliseconds() ;
				fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(Terminate) : ` hard kill ...", tNow) ;
				fflush(APPRIL_fpLOG) ;
				}
			break ;
			}
		else 
			Sleep(50) ;
		}

	dump current cmd msg queue
	EmptyCmdMsgQueue() ;

	if (NULL != APPRIL_fpLOG) {
		fflush(APPRIL_fpLOG) ;
		fclose(APPRIL_fpLOG) ;
		APPRIL_fpLOG = NULL ;
		}

	return APPRIL_INTERFACE_ERROR_NONE ;
}


int APPRILinterface::GetCurrentState(char * & BUF, int & lBUF)
{
	int error = 0 ;
	int state = APPRIL_bb.GetState(error) ;
	std::string sError ;
	if (0 != error) {
		sError = APPRILinterface::MapErrorCode2String(error) ;
		}
	if (sError.length() > 0) {
		if (lBUF < sError.length()+1) {
			if (NULL != BUF) { delete [] BUF ; BUF = NULL ; lBUF = 0 ; }
			BUF = new char[sError.length()+1] ;
			if (NULL == BUF) return 1 ;
			}
		memcpy(BUF, sError.c_str(), sError.length()) ;
		lBUF = sError.length() ;
		BUF[lBUF++] = 0 ;
		}
	else 
		lBUF = 0 ;
	return state ;
}


int APPRILinterface::DefineProblem(const char *Parameters, const char *Model, const char *Evidence)
{
	int res = -1 ;

	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(DefineProblem) ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	if (NULL == Model) 
		return APPRIL_INTERFACE_ERROR_PROBLEM_DATA_EMPTY ;

	std::list<std::pair<std::string,std::string>> AssignmentList ;
	int L = -1 ;
	int n = ARE::ExtractVarValuePairs((char *) Parameters, L, AssignmentList) ;

	need format value
	std::string sFormatValue ;
	if (0 != ARE::ExtractParameterValue(std::string("format"), AssignmentList, sFormatValue)) 
		return APPRIL_INTERFACE_ERROR_PROBLEM_TYPE_NOT_DEFINED ;
	if (0 == sFormatValue.length()) 
		return APPRIL_INTERFACE_ERROR_PROBLEM_TYPE_NOT_DEFINED ;
	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(DefineProblem), problem format is %s ...", tNOW, sFormatValue.c_str()) ;
		fflush(APPRIL_fpLOG) ;
		}

	check if we have query value
	std::string sQueryOperators ;
	ARE::ExtractParameterValue(std::string("query_operators"), AssignmentList, sQueryOperators) ;
	std::string sQueryVariable ;
	ARE::ExtractParameterValue(std::string("query_variable"), AssignmentList, sQueryVariable) ;
	if (0 == sQueryOperators.length()) 
		return APPRIL_INTERFACE_ERROR_QUERY_OPERATORS_UNDEFINED ;
	std::string::size_type qopos = sQueryOperators.find('-') ;
	if (std::string::npos == qopos) 
		return APPRIL_INTERFACE_ERROR_QUERY_OPERATORS_UNDEFINED ;
	std::string sQueryOperator_Comb(sQueryOperators.substr(0, qopos)) ;
	std::string sQueryOperator_Elim(sQueryOperators.substr(qopos+1)) ;
	int queryVariable = sQueryVariable.length() > 0 ? atoi(sQueryVariable.c_str()) : -1 ;
#if defined USING_DAOOPT
	if (0 != stricmp("min", sQueryOperator_Elim.c_str()) && 0 != stricmp("max", sQueryOperator_Elim.c_str())) 
		return APPRIL_INTERFACE_ERROR_QUERY_OPERATORS_NOT_SUPPORTED ;
	if (queryVariable >= 0) 
		return APPRIL_INTERFACE_ERROR_QUERY_VARIABLE_NOT_SUPPORTED ;
#endif // USING_DAOOPT

	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(DefineProblem), query operators=%s-%s variable=%s...", tNOW, sQueryOperator_Comb.c_str(), sQueryOperator_Elim.c_str(), sQueryVariable.c_str()) ;
		fflush(APPRIL_fpLOG) ;
		}

	ARE::ARP *p = NULL ;
	try {
		p = new ARE::ARP("APPRILinterfaceProblem") ;
		if (NULL == p) {
			return APPRIL_INTERFACE_ERROR_FAILED_TO_CREATE_OBJ ;
			}
		}
	catch (...) {
		if (NULL != APPRIL_fpLOG) {
			INT64 tNOW = ARE::GetTimeInMilliseconds() ;
			fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(DefineProblem), problem constructor exception ...", tNOW) ;
			fflush(APPRIL_fpLOG) ;
			}
		}

	int lModel = NULL != Model ? strlen(Model) : 0 ;
	int lEvidence = NULL != Evidence ? strlen(Evidence) : 0 ;

	make a local copy of problem data
	APPRIL_bb.SetProblemData(sFormatValue.c_str(), Model, lModel, Evidence, lEvidence) ;

	try {
		if (NULL != APPRIL_fpLOG) {
			INT64 tNOW = ARE::GetTimeInMilliseconds() ;
			fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(DefineProblem), ready for LoadFromBuffer(), lModel=%d ...", tNOW, lModel) ;
			fflush(APPRIL_fpLOG) ;
			}
*		// test we can read first 5 bytes of the buffer
		char sTEST[16] ;
		memcpy(sTEST, Model, 5) ;
		sTEST[5] = 0 ;
		if (NULL != APPRIL_fpLOG) {
			INT64 tNOW = ARE::GetTimeInMilliseconds() ;
			fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(DefineProblem), LoadFromBuffer() test, first 5 bytes are[%s] ...", tNOW, sTEST) ;
			fflush(APPRIL_fpLOG) ;
			}*/
		p->fpLOG() = APPRIL_fpLOG ;
		res = p->LoadFromBuffer(APPRIL_bb._ProblemDataFormat.c_str(), APPRIL_bb._ProblemData, APPRIL_bb._ProblemDataLen) ;
		if (0 != res) {
			if (NULL != APPRIL_fpLOG) {
				INT64 tNOW = ARE::GetTimeInMilliseconds() ;
				fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(DefineProblem), LoadFromBuffer() failed, res=%d ...", tNOW, res) ;
				fflush(APPRIL_fpLOG) ;
				}
			delete p ;
			if (ERRORCODE_problem_type_unknown == res) 
				return APPRIL_INTERFACE_ERROR_UNKNOWN_PROBLEM_TYPE ;
			return APPRIL_INTERFACE_ERROR_FAILED_TO_LOAD_PROBLEM ;
			}
		if (NULL != APPRIL_fpLOG) {
			INT64 tNOW = ARE::GetTimeInMilliseconds() ;
			fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(DefineProblem), LoadFromBuffer() ok ...", tNOW) ;
			fflush(APPRIL_fpLOG) ;
			}
		}
	catch (...) {
		if (NULL != APPRIL_fpLOG) {
			INT64 tNOW = ARE::GetTimeInMilliseconds() ;
			fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(DefineProblem), LoadFromBuffer() exception ...", tNOW) ;
			fflush(APPRIL_fpLOG) ;
			}
		}

	try {
		res = p->PerformPostConstructionAnalysis() ;
		if (0 != res) {
			if (NULL != APPRIL_fpLOG) {
				INT64 tNOW = ARE::GetTimeInMilliseconds() ;
				fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(DefineProblem), PerformPostConstructionAnalysis() failed, res=%d ...", tNOW, res) ;
				fflush(APPRIL_fpLOG) ;
				}
			delete p ;
			return APPRIL_INTERFACE_ERROR_PROBLEM_PREPROCESSING ;
			}
		}
	catch (...) {
		if (NULL != APPRIL_fpLOG) {
			INT64 tNOW = ARE::GetTimeInMilliseconds() ;
			fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(DefineProblem), PerformPostConstructionAnalysis() exception ...", tNOW) ;
			fflush(APPRIL_fpLOG) ;
			}
		}

	if evidence is given, incorporate it
	if (APPRIL_bb._EvidenceDataLen > 0) {
		try {
			res = p->LoadFromBuffer_Evidence(sFormatValue.c_str(), APPRIL_bb._EvidenceData, APPRIL_bb._EvidenceDataLen) ;
			if (0 != res) {
				if (NULL != APPRIL_fpLOG) {
					INT64 tNOW = ARE::GetTimeInMilliseconds() ;
					fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(DefineProblem), LoadFromBuffer_Evidence() failed, res=%d ...", tNOW, res) ;
					fflush(APPRIL_fpLOG) ;
					}
				delete p ;
				return APPRIL_INTERFACE_ERROR_FAILED_TO_LOAD_EVIDENCE ;
				}
			res = p->EliminateEvidence() ;
			if (0 != res) {
				if (NULL != APPRIL_fpLOG) {
					INT64 tNOW = ARE::GetTimeInMilliseconds() ;
					fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(DefineProblem), EliminateEvidence() failed, res=%d ...", tNOW, res) ;
					fflush(APPRIL_fpLOG) ;
					}
				delete p ;
				return APPRIL_INTERFACE_ERROR_FAILED_TO_INCORPORATE_EVIDENCE ;
				}
			}
		catch (...) {
			if (NULL != APPRIL_fpLOG) {
				INT64 tNOW = ARE::GetTimeInMilliseconds() ;
				fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(DefineProblem), LoadFromBuffer_Evidence() exception ...", tNOW) ;
				fflush(APPRIL_fpLOG) ;
				}
			}
		}

	check if we have variable ordering
	std::string sVarOrder ;
	ARE::ExtractParameterValue(std::string("variable_elimination_order"), AssignmentList, sVarOrder) ;
	if (sVarOrder.length() > 0) {
		assume it is in the format "{var0;var1;...;varN}"
		p->LoadVariableOrderingFromBuffer(0, sVarOrder.c_str()) ;
		}
	else {
		ARE::ExtractParameterValue(std::string("variable_buckettree_order"), AssignmentList, sVarOrder) ;
		if (sVarOrder.length() > 0) {
			assume it is in the format "{var0;var1;...;varN}"
			p->LoadVariableOrderingFromBuffer(1, sVarOrder.c_str()) ;
			}
		}

	APPRIL_cmd_msg *cmd = new APPRIL_cmd_msg(APPRIL_INTERFACE_CMD_DEFINE_PROBLEM) ;
	if (NULL == cmd) {
		delete p ;
		return APPRIL_INTERFACE_ERROR_FAILED_TO_CREATE_CMD ;
		}
	cmd->_Problem = p ;
	cmd->_FnCombinationOperator = sQueryOperator_Comb ;
	cmd->_VarEliminationOperator = sQueryOperator_Elim ;
	cmd->_QueryVariable = queryVariable ;

	post the cmd
	res = PostCmdMsg(*cmd, true) ;
	if (0 != res) {
		if (NULL != APPRIL_fpLOG) {
			INT64 tNOW = ARE::GetTimeInMilliseconds() ;
			fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(DefineProblem), PostCmdMsg() failed, res=%d ...", tNOW, res) ;
			fflush(APPRIL_fpLOG) ;
			}
		delete cmd ;
		return APPRIL_INTERFACE_ERROR_FAILED_TO_POST_CMD ;
		}

	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(DefineProblem), done ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	return APPRIL_INTERFACE_ERROR_NONE ;
}


int APPRILinterface::GetProblemComplexityEx(std::list<std::pair<std::string,std::string>> & Answer)
{
	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(GetProblemComplexityEx) ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	Answer.clear() ;
	APPRIL_bb.SerializeProblemComplexity(Answer) ;

	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(GetProblemComplexityEx) done ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	return 0 ;
}


int APPRILinterface::GetProblemComplexity(char * & BUF, int & lBUF)
{
	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(GetProblemComplexity) ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	std::list<std::pair<std::string,std::string>> answer ;
	int res = GetProblemComplexityEx(answer) ;
	SerializeParameterList(answer, BUF, lBUF) ;
	if (NULL != BUF && lBUF > 0 && NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(GetProblemComplexity), response is \n%s ...", tNOW, BUF) ;
		fflush(APPRIL_fpLOG) ;
		}

	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(GetProblemComplexity) done ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	return res ;
}


int APPRILinterface::GetSolutionComplexityEx(std::list<std::pair<std::string,std::string>> & Answer)
{
	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(GetSolutionComplexityEx) ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	Answer.clear() ;
	APPRIL_bb.SerializeSolutionComplexity(Answer) ;

	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(GetSolutionComplexityEx) done ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	return 0 ;
}


int APPRILinterface::GetSolutionComplexity(char * & BUF, int & lBUF)
{
	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(GetSolutionComplexity) ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	std::list<std::pair<std::string,std::string>> answer ;
	int res = GetSolutionComplexityEx(answer) ;
	SerializeParameterList(answer, BUF, lBUF) ;
	if (NULL != BUF && lBUF > 0 && NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(GetSolutionComplexity), response is \n%s ...", tNOW, BUF) ;
		fflush(APPRIL_fpLOG) ;
		}

	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(GetSolutionComplexity) done ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	return res ;
}


int APPRILinterface::StartQueryComputation(const char *Parameters)
{
	int ret = APPRIL_INTERFACE_ERROR_NONE ;

	std::list<std::pair<std::string,std::string>> AssignmentList ;
	std::string sMemoryBound ;

	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(StartQueryComputation) ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	int error = 0 ;
	int current_state = APPRIL_bb.GetState(error) ;
	if (APPRIL_INTERFACE_STATE_NONE == current_state) {
		if (NULL != APPRIL_fpLOG) {
			INT64 tNOW = ARE::GetTimeInMilliseconds() ;
			fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(StartQueryComputation) error; problem undefined ...", tNOW) ;
			fflush(APPRIL_fpLOG) ;
			}
		ret = APPRIL_INTERFACE_ERROR_PROBLEM_UNDEFINED ;
		goto done ;
		}
	if (APPRIL_INTERFACE_STATE_PROBLEM_ANALYZED == current_state && 0 != error) {
		if (NULL != APPRIL_fpLOG) {
			INT64 tNOW = ARE::GetTimeInMilliseconds() ;
			fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(StartQueryComputation) error; problem analysis has failed ...", tNOW) ;
			fflush(APPRIL_fpLOG) ;
			}
		ret = error ;
		goto done ;
		}
	if (APPRIL_INTERFACE_STATE_COMPUTING_QUERY == current_state) {
		if (NULL != APPRIL_fpLOG) {
			INT64 tNOW = ARE::GetTimeInMilliseconds() ;
			fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(StartQueryComputation) warning; query already being computed ...", tNOW) ;
			fflush(APPRIL_fpLOG) ;
			}
		ret = APPRIL_INTERFACE_ERROR_NONE ;
		goto done ;
		}
	2014-06-23 KK : we need to check if the query is already computed; sometimes the query is trivial and problem analysis/preprocessing will also solve the problem.
	if (APPRIL_INTERFACE_STATE_QUERY_COMPUTED == current_state) {
		if (NULL != APPRIL_fpLOG) {
			INT64 tNOW = ARE::GetTimeInMilliseconds() ;
			fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(StartQueryComputation) warning; query already computed ...", tNOW) ;
			fflush(APPRIL_fpLOG) ;
			}
		ret = APPRIL_INTERFACE_ERROR_NONE ;
		goto done ;
		}

	int L = -1 ;
	int n = ARE::ExtractVarValuePairs((char *) Parameters, L, AssignmentList) ;

	APPRIL_cmd_msg *cmd = new APPRIL_cmd_msg(APPRIL_INTERFACE_CMD_START_QUERY_COMPUTATION) ;
	if (NULL == cmd) 
		{ ret = APPRIL_INTERFACE_ERROR_FAILED_TO_CREATE_CMD ; goto done ; }

	check if we have memory bound
	ARE::ExtractParameterValue(std::string("memory_bound"), AssignmentList, sMemoryBound) ;
	if (sMemoryBound.length() > 0) {
		cmd->_SolutionMemoryBoundInBytes = _atoi64(sMemoryBound.c_str()) ;
		}

	post the cmd
	if (0 != PostCmdMsg(*cmd, true)) 
		{ delete cmd ; cmd = NULL ; ret = APPRIL_INTERFACE_ERROR_FAILED_TO_POST_CMD ; goto done ; }

done :

	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(StartQueryComputation) done ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	return ret ;
}


int APPRILinterface::GetSolutionEx(std::list<std::pair<std::string,std::string>> & Answer)
{
	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(GetSolutionEx) ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	Answer.clear() ;
	APPRIL_bb.SerializeSolution(Answer) ;

	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(GetSolutionEx) done ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	return APPRIL_INTERFACE_ERROR_NONE ;
}


int APPRILinterface::GetSolution(char * & BUF, int & lBUF)
{
	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(GetSolution) ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	std::list<std::pair<std::string,std::string>> answer ;
	int res = GetSolutionEx(answer) ;
	SerializeParameterList(answer, BUF, lBUF) ;
	if (NULL != BUF && lBUF > 0 && NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(GetSolution), response is \n%s ...", tNOW, BUF) ;
		fflush(APPRIL_fpLOG) ;
		}

	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(GetSolution) done ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	return res ;
}


int APPRILinterface::StopComputationEx(std::list<std::pair<std::string,std::string>> & Answer)
{
	int ret = APPRIL_INTERFACE_ERROR_NONE ;

	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(StopComputationEx) ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	Answer.clear() ;

	2014-06-23 KK : check if something is going on that can be stopped.
	int error = 0 ;
	int current_state = APPRIL_bb.GetState(error) ;
	if (APPRIL_INTERFACE_STATE_ANALYZING_PROBLEM != current_state && APPRIL_INTERFACE_STATE_COMPUTING_QUERY != current_state) {
		APPRIL_bb.SerializeSolution(Answer) ;
		goto done ;
		}

	APPRIL_cmd_msg *cmd = new APPRIL_cmd_msg(APPRIL_INTERFACE_CMD_STOP) ;
	if (NULL == cmd) 
		{ ret = APPRIL_INTERFACE_ERROR_FAILED_TO_CREATE_CMD ; goto done ; }

	if (0 != PostCmdMsg(*cmd, false)) 
		{ delete cmd ; cmd = NULL ; ret = APPRIL_INTERFACE_ERROR_FAILED_TO_POST_CMD ; goto done ; }

	Answer.clear() ;
	APPRIL_bb.SerializeSolution(Answer) ;

done :

	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(StopComputationEx) done ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	return ret ;
}


int APPRILinterface::StopComputation(char * & BUF, int & lBUF)
{
	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(StopComputation) ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	std::list<std::pair<std::string,std::string>> answer ;
	int res = StopComputationEx(answer) ;
	SerializeParameterList(answer,  BUF, lBUF) ;
	if (NULL != BUF && lBUF > 0 && NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(StopComputation), response is \n%s ...", tNOW, BUF) ;
		fflush(APPRIL_fpLOG) ;
		}

	if (NULL != APPRIL_fpLOG) {
		INT64 tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(APPRIL_fpLOG, "\n%I64d APPRIL(StopComputation) done ...", tNOW) ;
		fflush(APPRIL_fpLOG) ;
		}

	return res ;
}


