#ifndef RL_DIST_DIVERSITY_DEPENDENT_FBD_H
#define RL_DIST_DIVERSITY_DEPENDENT_FBD_H
#include "DiversityDependentFossilizedBirthDeathProcess.h"
#include "RlTimeTree.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
class Dist_diversityDependentFBD : public TypedDistribution<TimeTree> {
public:
    Dist_diversityDependentFBD();
    Dist_diversityDependentFBD* clone(void) const override;
    static const std::string& getClassType(void);
    static const TypeSpec& getClassTypeSpec(void);
    const TypeSpec& getTypeSpec(void) const override;
    std::string getDistributionFunctionName(void) const override;
    std::vector<std::string> getDistributionFunctionAliases(void) const override;
    const MemberRules& getParameterRules(void) const override;
    RevBayesCore::DiversityDependentFossilizedBirthDeathProcess* createDistribution(void) const override;
protected:
    void setConstParameter(const std::string& name, const RevPtr<const RevVariable>& var) override;
private:
    RevPtr<const RevVariable> origin, lambda0, alpha, mu, psi, rho;
    RevPtr<const RevVariable> occurrenceSpecies, occurrenceAges, sampledExtant;
    RevPtr<const RevVariable> initialTree, condition, maxHidden, tolerance;
};
}
#endif
