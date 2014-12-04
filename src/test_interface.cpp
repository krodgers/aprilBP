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
  // set up args so don't have to worry about command line parameters for debugging
  // argc = 6;
  // char** params = (char**) malloc(sizeof(char*) * 5);
  // params[0] = "linux/build/bp_interface"; 
  // params[1] = "-f";
  // params[2] = "../data/test.uai";
  // params[3] = "-T";
  // params[4] = "PR";
  // params[5] = "--verbose";
  
  
  printf("Test Iterface\n");
  printf("Starting...\n");
  BpInterface bpi;
  printf("Initializing\n");
  bool success = bpi.initialize(argc,argv);
 //free(params);
  if(!success){
    printf("Failed Initializing\n");
    return 0;
  }
  
  printf("Analyzing\n");
  int time, memory;
  bpi.estimateComplexity(time, memory);
  printf("Estimated %d time and %d memory\n", time, memory);

  printf("Beginning Inference\n");
  bpi.runInference(120);
  printf("Getting solution\n");
  double ans =  bpi.getPRSolution();
  printf("Solution returned: %.4f\n",ans);
  printf("Quitting\n");
  return 0;
}

void writePR(const char* outfile, double logZ) {
  ofstream os(outfile);
  //os.precision(8); os.setf(ios::fixed,ios::floatfield);
  //os<<"PR\n1\n"<<logZ/c_log10<<"\n";		// !!! 2011 PIC version: need "1" evidence
  os<<"PR\n"<<logZ/c_log10<<"\n";
  os.close();
  std::cout<<"Wrote PR : "<<logZ/c_log10<<"\n";
}

void writeMAR(const char* outfile, mex::vector<Factor>& fs) {
  ofstream os(outfile);
  //os<<"MAR\n1\n";		// !!! 2011 PIC version: need "1" evidence
  os<<"MAR\n";
  os<<fs.size()<<" ";
  for (size_t f=0;f<fs.size();++f) {
    os<<fs[f].nrStates()<<" ";
    for (size_t i=0;i<fs[f].nrStates();++i) os<<fs[f][i]<<" ";
  }
  os<<"\n";
  os.close();
  std::cout<<"Wrote MAR\n";
}

