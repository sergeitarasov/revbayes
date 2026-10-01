#ifndef RB_DIVERSITY_DEPENDENT_FOSSILIZED_BIRTH_DEATH_PROCESS_H
#define RB_DIVERSITY_DEPENDENT_FOSSILIZED_BIRTH_DEATH_PROCESS_H

#include "DiversityDependentFbdLikelihood.h"
#include "RbVector.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include <map>
#include <set>
#include <string>
#include <vector>

namespace RevBayesCore {
/** Joint density of an oriented extended budding tree and assigned fossil records.
 * Child zero retains ancestral identity. Tips are true death/present endpoints;
 * characters must attach through fnSpeciesObservationTree, not to these tips.
 * Living/extinct status is fixed by initialTree. Hidden species are marginalized.
 */
class DiversityDependentFossilizedBirthDeathProcess : public TypedDistribution<Tree> {
public:
    DiversityDependentFossilizedBirthDeathProcess(
        const TypedDagNode<double>* origin, const TypedDagNode<double>* lambda0,
        const TypedDagNode<double>* alpha, const TypedDagNode<double>* mu,
        const TypedDagNode<double>* psi, const TypedDagNode<double>* rho,
        const std::vector<std::string>& occurrenceSpecies,
        const TypedDagNode<RbVector<double> >* occurrenceAges,
        const std::vector<std::string>& sampledExtant, const Tree& initialTree,
        const std::string& condition, const DDFBD::Settings& settings);
    DiversityDependentFossilizedBirthDeathProcess* clone(void) const override;
    double computeLnProbability(void) override;
    void redrawValue(void) override;
    void redrawValue(SimulationCondition condition) override;
protected:
    void swapParameterInternal(const DagNode* oldP, const DagNode* newP) override;
private:
    const TypedDagNode<double>* originAge;
    const TypedDagNode<double>* birthRate;
    const TypedDagNode<double>* feedback;
    const TypedDagNode<double>* deathRate;
    const TypedDagNode<double>* fossilRate;
    const TypedDagNode<double>* extantSampling;
    const TypedDagNode<RbVector<double> >* fossilAges;
    std::vector<std::string> fossilSpecies;
    std::set<std::string> sampledLiving;
    std::map<std::string, bool> livingStatus;
    Tree initial;
    std::string conditioning;
    DDFBD::Settings numerical;
};
}
#endif
