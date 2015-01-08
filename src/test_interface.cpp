#include <stdio.h>
#include "BpInterface.h"
#include <fstream>
#include <iostream>

/*

  A testing file for testing interactions with BpInterface
 */


using namespace std;
using lgbp::BpInterface;

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

int main(int argc, char** argv){
  
  printf("Test Iterface\n");
  printf("Starting...\n");
  BpInterface bpi;
  printf("Initializing\n");
  char t[] = "PR";
  char prob[] = {"../data/test.uai"};
  char order[] = {"../data/eliminationOrder"};
  //  bool success = bpi.initialize(1800, t, prob, order, NULL, true);
  bool success = bpi.initialize(argc, argv);

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
  bpi.runInference();
  printf("Getting solution\n");
  double ans;
  mex::vector<Factor> MarSltn(9);
  bpi.getSolution(MarSltn);
  //  bpi.getSolution(ans);
  //  printf("Solution returned:%g \n",ans);
  // std::cout<<ans<<"\n";
  writeMAR("../results/Interface_MAR", MarSltn);
  printf("Quitting\n");
  return 0;
}


