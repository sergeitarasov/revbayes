#ifndef RL_DIST_DD_SPECIMEN_FBD_H
#define RL_DIST_DD_SPECIMEN_FBD_H
#include "DiversityDependentSpecimenBirthDeathProcess.h"
#include "RlTimeTree.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
class Dist_DDFBDP : public TypedDistribution<TimeTree> {
public:
    Dist_DDFBDP();
    Dist_DDFBDP* clone(void) const override;
    static const std::string& getClassType(void);
    static const TypeSpec& getClassTypeSpec(void);
    const TypeSpec& getTypeSpec(void) const override;
    std::string getDistributionFunctionName(void) const override;
    std::vector<std::string> getDistributionFunctionAliases(void) const override;
    const MemberRules& getParameterRules(void) const override;
    RevBayesCore::DiversityDependentSpecimenBirthDeathProcess* createDistribution(void) const override;
protected:
    void setConstParameter(const std::string& name, const RevPtr<const RevVariable>& var) override;
private:
    RevPtr<const RevVariable> origin, lambda0, alpha, mu, psi, rho;
    RevPtr<const RevVariable> taxa, agePrecision;
    RevPtr<const RevVariable> initialTree, condition, maxHidden, tolerance;
};
}
#endif
