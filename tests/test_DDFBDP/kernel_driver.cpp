#include "DiversityDependentFbdLikelihood.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <exception>
using namespace RevBayesCore::DDFBD;
int main(int argc,char** argv) {
 try {
  if (argc<10 || (argc-10)%2) return 2;
  Parameters p{std::stod(argv[1]),std::stod(argv[2]),std::stod(argv[3]),std::stod(argv[4]),std::stod(argv[5])};
  Settings s; s.maxHidden=std::stoul(argv[7]); s.tolerance=std::stod(argv[8]);
  std::vector<SpecimenEvent> e;
  for(int i=10;i<argc;i+=2) {
   std::string t=argv[i+1];
   e.push_back({std::stod(argv[i]),t=="B"?SpecimenEventType::Birth:t=="A"?SpecimenEventType::AncestorSample:SpecimenEventType::TerminalSample});
  }
  std::cout<<std::setprecision(17)<<specimenLogLikelihood(std::stod(argv[6]),e,std::stoul(argv[9]),p,s)<<'\n';
 } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 2; }
}
