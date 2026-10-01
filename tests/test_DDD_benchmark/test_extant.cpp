#include "DiversityDependentFbdLikelihood.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace RevBayesCore::DDFBD;
void equal(double a,double b) {
    if(!std::isfinite(a)||!std::isfinite(b)||std::abs(a-b)>1e-9)
        throw std::runtime_error("Extant/FBD reduction or pure-birth formula mismatch");
}
int main() {
    const std::vector<double> ages{3,2,1};
    Settings settings{128,1e-13};
    unsigned checks=0;
    for(BirthModel model : {BirthModel::Exponential,BirthModel::DDDLinear,BirthModel::DDDPower}) {
        Parameters p{.8,.2,.2,0,1,model,8};
        for(bool crown : {false,true}) {
            const double origin=crown?3:4;
            std::vector<Event> events{{3,EventType::Birth},{2,EventType::Birth},{1,EventType::Birth}};
            for(unsigned i=0;i<4;++i) events.push_back({0,EventType::Protect});
            const double fbd=logLikelihood(origin,events,0,4,0,p,settings);
            const double dd=extantTreeLogLikelihood(origin,ages,crown,0,false,p,settings);
            equal(dd,fbd-(crown?std::log(p.birth(1)):0)); ++checks;
            equal(extantTreeLogLikelihood(origin,ages,crown,0,true,p,settings),dd+std::lgamma(4.)); ++checks;
        }
    }
    for(double alpha : {0.,.2,2.}) for(bool crown : {false,true}) {
        Parameters p{.8,alpha,0,0,1};
        double expected=0;
        for(unsigned k=crown?2:1;k<=4;++k) {
            expected-=k*p.birth(k); // each interval lasts one unit
            if(k<4) expected+=std::log(p.birth(k));
        }
        equal(extantTreeLogLikelihood(crown?3:4,ages,crown,0,false,p,settings),expected); ++checks;
    }
    bool rejected=false;
    try { extantTreeLogLikelihood(2,ages,false,0,false,Parameters{.8,.2,.2,0,1},settings); }
    catch(const std::invalid_argument&) { rejected=true; }
    if(!rejected) throw std::runtime_error("Invalid stem age accepted");
    ++checks;
    std::cout << checks << " extant reduction/analytic/domain checks passed\n";
}
