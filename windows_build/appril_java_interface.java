

public class appril_java_interface{

    // load the .dll
    static{
	System.loadLibrary("LBP");
    }

    // callable Methods

    // return value = error string
    public native   char[] MapErrorCode2String(int ErrorCode) ;

    // return value = error code
    public native int Initialize(boolean logStuff) ;
    // return value = error code
    // defaults to logging
    public native int Initialize() ;


    // return value = error code
    public native int Terminate() ;

    // return value = current state
    public native int GetCurrentState(char[] BUF, int lBUF) ;

    // return value = error code
    public native int DefineProblem(  char[] Parameters,   char[] Model,   char[] Evidence) ;

    // return value = error code
    //    public native int GetProblemComplexityEx(std::list<std::pair<std::string,std::string> >  Answer) ;
    public native int GetProblemComplexity(char[] BUF, int lBUF) ;

    // return value = error code
    //public native int GetSolutionComplexityEx(std::list<std::pair<std::string,std::string> >  Answer) ;
    public native int GetSolutionComplexity(char[] BUF, int lBUF) ;

    // return value = error code
    public native int StartQueryComputation(  char[] Parameters) ;

    // return value = error code
    //public native int GetSolutionEx(std::list<std::pair<std::string,std::string> >  Answer) ;
    public native int GetSolution(char[] BUF, int  lBUF) ;

    // return value = error code
    //public native int StopComputationEx(std::list<std::pair<std::string,std::string> >  Answer) ;
    public native int StopComputation(char[] BUF, int  lBUF) ;


}
