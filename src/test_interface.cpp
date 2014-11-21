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
  bpi.initialize(argc,argv);
  printf("Analyzing\n");
  bpi.estimateComplexity();
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

