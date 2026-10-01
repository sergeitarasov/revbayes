// Compile-only regression: FBDSP must override the const virtual used by
// AbstractRootedTreeDistribution::computeLnProbability().
#include "FossilizedBirthDeathSpeciationProcess.h"
#include <type_traits>

struct FbdspDispatchProbe : RevBayesCore::FossilizedBirthDeathSpeciationProcess {
    using FossilizedBirthDeathSpeciationProcess::computeLnProbabilityDivergenceTimes;
};

static_assert(std::is_same_v<
    decltype(&FbdspDispatchProbe::computeLnProbabilityDivergenceTimes),
    double (RevBayesCore::FossilizedBirthDeathSpeciationProcess::*)() const>,
    "FBDSP must override the const divergence-times likelihood virtual");
