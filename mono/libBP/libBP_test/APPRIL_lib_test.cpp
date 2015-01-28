#include <stdio.h>
#include "APPRIL_INTERFACE.h"
#include <fstream>
#include <iostream>

/*

  A testing file for testing interactions with BpInterface
 */


using namespace std;




int main(int argc, char** argv){
  
  printf("Test Iterface\n");
  printf("Starting...\n");
  printf("Initializing\n");
  int success = APPRILinterface::Initialize(true);
 #ifdef MAR
 char t[] = "MAR";
  
 #else
 char t[] = "PR";
  
 #endif
 //char prob[] = {"/home/krodgers/Documents/Research/aprilBP/data/test.uai"};
 // char order[] = {"/home/krodgers/Documents/Research/aprilBP/data/data/eliminationOrder"};
 char prob[] = {"/home/krodgers/Documents/Research/aprilBP/data/MAR_tests/uai08_test3.uai"};


  if(!success){
    printf("Failed Initializing\n");
 	return 0;
  }
  
  printf("Analyzing\n");
  int time, memory;
  printf("Estimated %d time and %d memory\n", time, memory);

  printf("Beginning Inference\n");
  printf("Getting solution\n");

  #ifdef MAR
  mex::vector<Factor> MarSltn(9);
  #else
  double ans;
  #endif

  /// DELETE ME ////
  printf("Returned from get solution\n");
  //////////////////
  
  printf("Solution returned:%g \n",ans);
  printf("Quitting\n");

  int res;
  char* buf;
  res = APPRILinterface::GetSolution(buf, 0);
  if(!res)
    printf("Failed getting solution\n");
  return 0;
}



