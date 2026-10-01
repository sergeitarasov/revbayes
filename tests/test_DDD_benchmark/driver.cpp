#include "DiversityDependentFbdLikelihood.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
using namespace RevBayesCore::DDFBD;
int main(int argc,char** argv) {
    try {
        if(argc<11) throw std::invalid_argument("model lambda mu alpha K crown condition cap startAge branchingAges...");
        Parameters p{std::stod(argv[2]),std::stod(argv[4]),std::stod(argv[3]),0,1};
        const std::string model=argv[1];
        if(model=="DDDlinear") p.birthModel=BirthModel::DDDLinear;
        else if(model=="DDDpower") p.birthModel=BirthModel::DDDPower;
        else if(model!="exponential") throw std::invalid_argument("Invalid rate model");
        p.capacity=std::stod(argv[5]);
        Settings settings{std::stoull(argv[8]),1e-13};
        std::vector<double> ages;
        for(int i=10;i<argc;++i) ages.push_back(std::stod(argv[i]));
        std::cout << std::setprecision(17) << extantTreeLogLikelihood(std::stod(argv[9]),ages,
            std::stoi(argv[6])!=0,std::stoul(argv[7]),false,p,settings) << '\n';
    } catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
