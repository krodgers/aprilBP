#include <stdio.h>
#include "BPinterface.h"
#include <fstream>
#include <iostream>

/*

  A testing file for testing interactions with BpInterface
 */


using namespace std;
using lgbp::BpInterface;




int main(int argc, char** argv){
  
  printf("Test Iterface\n");
  printf("Starting...\n");
  BpInterface bpi;
  printf("Initializing\n");
 #ifdef MAR
 char t[] = "MAR";
  
 #else
 char t[] = "PR";
  
 #endif
 //char prob[] = {"/home/krodgers/Documents/Research/aprilBP/data/test.uai"};
 // char order[] = {"/home/krodgers/Documents/Research/aprilBP/data/data/eliminationOrder"};
 char prob[] = {"/home/krodgers/Documents/Research/aprilBP/data/MAR_tests/uai08_test3.uai"};
bool success = bpi.initialize(1800, t, prob, NULL, NULL, true);
 // bool success = bpi.initialize(1800, t, prob, order, NULL, true);
  // bool success = bpi.initialize(argc, argv);

 
  if(!success){
    printf("Failed Initializing\n");
 	return 0;
  }
  
  printf("Analyzing\n");
  int time, memory;
  bpi.estimateComplexity(time, memory);
  printf("Estimated %d time and %d memory\n", time, memory);

  printf("Beginning Inference\n");
  bpi.runInference();
  printf("Getting solution\n");

  #ifdef MAR
  mex::vector<Factor> MarSltn(9);
  bool res = bpi.getSolution(MarSltn); // gets solution and writes to file
  #else
  double ans;
  bool res = bpi.getSolution(ans);
  #endif

  /// DELETE ME ////
  printf("Returned from get solution\n");
  //////////////////
  
  printf("Solution returned:%g \n",ans);
  printf("Quitting\n");

  // int dummyVar;
  //std::cin >> dummyVar;
  if(!res)
    printf("Failed getting solution\n");
  return 0;
}


