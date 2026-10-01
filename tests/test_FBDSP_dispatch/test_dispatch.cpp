// Runtime dispatch regression. Override the range evaluator with a sentinel so
// this tests virtual wiring independently of the range model's mathematics.
#include "ConstantNode.h"
#include "FossilizedBirthDeathSpeciationProcess.h"
#include "TimeInterval.h"
#include "RbMathCombinatorialFunctions.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

class DispatchProbe : public RevBayesCore::FossilizedBirthDeathSpeciationProcess {
public:
    using FossilizedBirthDeathSpeciationProcess::FossilizedBirthDeathSpeciationProcess;
    bool evaluated_ranges = false;
protected:
    double computeLnProbabilityRanges(bool = false) override {
        evaluated_ranges = true;
        return 13.0;
    }
    double computeLnProbabilityTimes() const override { return 17.0; }
};

int main() {
    using namespace RevBayesCore;
    auto constant = [](const char* name, double x) {
        return new ConstantNode<double>(name, new double(x));
    };
    Taxon a("a"), b("b");
    a.addOccurrence(TimeInterval(0.0, 0.0));
    b.addOccurrence(TimeInterval(0.0, 0.0));
    DispatchProbe distribution(constant("origin", 3.0), constant("lambda", 0.4),
        constant("mu", 0.1), constant("psi", 0.2), constant("rho", 1.0),
        constant("lambda_a", 0.0), constant("beta", 0.0), nullptr,
        "time", std::vector<Taxon>{a, b}, true, false);
    AbstractRootedTreeDistribution& base = distribution;
    const double actual = base.computeLnProbability();
    const double expected = 30.0 - RbMath::lnFactorial(2);
    if (!distribution.evaluated_ranges || std::abs(actual - expected) > 1e-12) {
        std::cerr << "FBDSP range dispatch failed: " << actual << " expected " << expected << '\n';
        return 1;
    }
    std::cout << "PASS: public tree likelihood dispatch includes FBDSP ranges.\n";
}
