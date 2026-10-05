#ifndef RB_DD_SPECIMEN_BIRTH_DEATH_PROCESS_H
#define RB_DD_SPECIMEN_BIRTH_DEATH_PROCESS_H
#include "AbstractBirthDeathProcess.h"
#include "DiversityDependentFbdLikelihood.h"
namespace RevBayesCore {
class DiversityDependentSpecimenBirthDeathProcess : public AbstractBirthDeathProcess {
public:
    DiversityDependentSpecimenBirthDeathProcess(const TypedDagNode<double>* origin,
        const TypedDagNode<double>* lambda, const TypedDagNode<double>* alpha,
        const TypedDagNode<double>* mu, const TypedDagNode<double>* psi,
        const TypedDagNode<double>* rho, const std::vector<Taxon>& taxa,
        Tree* initial, const std::string& condition, const DDFBD::Settings& settings,
        std::int64_t agePrecision);
    DiversityDependentSpecimenBirthDeathProcess* clone() const override;
    bool allowsSA() override { return true; }
    void redrawValue(SimulationCondition c) override;
    void redrawValue() override;
    void setValue(Tree* tree, bool force=false) override;
protected:
    double computeLnProbabilityDivergenceTimes() const override;
    double computeLnProbabilityTimes() const override;
    double lnProbTreeShape() const override;
    double lnProbNumTaxa(size_t,double,double,bool) const override;
    double pSurvival(double,double) const override;
    double simulateDivergenceTime(double,double) const override;
    std::vector<double> simulateDivergenceTimes(size_t,double,double,double,bool) const override;
    void swapParameterInternal(const DagNode*,const DagNode*) override;
private:
    const TypedDagNode<double> *lambda, *alpha, *mu, *psi, *rho;
    DDFBD::Settings settings;
    std::int64_t precision;
};
}
#endif
