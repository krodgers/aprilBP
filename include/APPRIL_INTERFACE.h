#ifndef APPRIL_INTERFACE_HXX_INCLUDED
#define APPRIL_INTERFACE_HXX_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <list>

#include <APPRIL_Export.hxx>

#define APPRIL_INTERFACE_VERSION_STRING "1.0.6"
// 1.0.6 - 2014-07-01 KK : first released version with VE/DAOOPT

#define APPRIL_INTERFACE_ERROR_NONE										0
#define APPRIL_INTERFACE_ERROR_UNKNOWN									1
#define APPRIL_INTERFACE_ERROR_SYSTEM									2
#define APPRIL_INTERFACE_ERROR_NOT_IMPLEMENTED							3
#define APPRIL_INTERFACE_ERROR_MEMORY_ALLOCATION						4
#define APPRIL_INTERFACE_ERROR_TIMEOUT									5
#define APPRIL_INTERFACE_ERROR_FAILED_TO_CREATE_OBJ						6
#define APPRIL_INTERFACE_ERROR_FAILED_TO_CREATE_CMD						7
#define APPRIL_INTERFACE_ERROR_FAILED_TO_POST_CMD						8
#define APPRIL_INTERFACE_ERROR_PROBLEM_TYPE_NOT_DEFINED					9
#define APPRIL_INTERFACE_ERROR_PROBLEM_DATA_EMPTY						10
#define APPRIL_INTERFACE_ERROR_UNKNOWN_PROBLEM_TYPE						11
#define APPRIL_INTERFACE_ERROR_PROBLEM_UNDEFINED						12
#define APPRIL_INTERFACE_ERROR_FAILED_TO_LOAD_PROBLEM					13
#define APPRIL_INTERFACE_ERROR_PROBLEM_PREPROCESSING					14
#define APPRIL_INTERFACE_ERROR_FAILED_TO_START_CVO						15
#define APPRIL_INTERFACE_ERROR_FAILED_TO_ANALYSE_PROBLEM				16
#define APPRIL_INTERFACE_ERROR_QUERY_UNDEFINED							17
#define APPRIL_INTERFACE_ERROR_QUERY_OPERATORS_UNDEFINED				18
#define APPRIL_INTERFACE_ERROR_QUERY_OPERATORS_NOT_SUPPORTED			19
#define APPRIL_INTERFACE_ERROR_QUERY_VARIABLE_NOT_SUPPORTED				20
#define APPRIL_INTERFACE_ERROR_QUERY_UNKNOWN							21
#define APPRIL_INTERFACE_ERROR_FAILED_TO_CONSTRUCT_SOLUTION_WORKSPACE	22
#define APPRIL_INTERFACE_ERROR_MEMORY_BOUND_EXCEEDED					23
#define APPRIL_INTERFACE_ERROR_FAILED_TO_LAUNCH_THREAD					24
#define APPRIL_INTERFACE_ERROR_FAILED_TO_LOAD_EVIDENCE					25
#define APPRIL_INTERFACE_ERROR_FAILED_TO_INCORPORATE_EVIDENCE			26

#define APPRIL_INTERFACE_STATE_NONE								0
#define APPRIL_INTERFACE_STATE_ANALYZING_PROBLEM				1
#define APPRIL_INTERFACE_STATE_PROBLEM_ANALYZED					2
#define APPRIL_INTERFACE_STATE_COMPUTING_QUERY					3
#define APPRIL_INTERFACE_STATE_QUERY_COMPUTED					4
#define APPRIL_INTERFACE_STATE_EXIT_COMPLETED					99

namespace APPRILinterface
{

// return value = error string
APPRIL_EXP_FUNC const char *MapErrorCode2String(int ErrorCode) ;

// return value = error code
APPRIL_EXP_FUNC int Initialize(bool logStuff = true) ;

// return value = error code
APPRIL_EXP_FUNC int Terminate(void) ;

// return value = current state
APPRIL_EXP_FUNC int GetCurrentState(char * & BUF, int & lBUF) ;

// return value = error code
APPRIL_EXP_FUNC int DefineProblem(const char *Parameters, const char *Model, const char *Evidence) ;

// return value = error code
APPRIL_EXP_FUNC int GetProblemComplexityEx(std::list<std::pair<std::string,std::string>> & Answer) ;
APPRIL_EXP_FUNC int GetProblemComplexity(char * & BUF, int & lBUF) ;

// return value = error code
APPRIL_EXP_FUNC int GetSolutionComplexityEx(std::list<std::pair<std::string,std::string>> & Answer) ;
APPRIL_EXP_FUNC int GetSolutionComplexity(char * & BUF, int & lBUF) ;

// return value = error code
APPRIL_EXP_FUNC int StartQueryComputation(const char *Parameters) ;

// return value = error code
APPRIL_EXP_FUNC int GetSolutionEx(std::list<std::pair<std::string,std::string>> & Answer) ;
APPRIL_EXP_FUNC int GetSolution(char * & BUF, int & lBUF) ;

// return value = error code
APPRIL_EXP_FUNC int StopComputationEx(std::list<std::pair<std::string,std::string>> & Answer) ;
APPRIL_EXP_FUNC int StopComputation(char * & BUF, int & lBUF) ;

} // namespace APPRILinterface

#endif // APPRIL_INTERFACE_HXX_INCLUDED
