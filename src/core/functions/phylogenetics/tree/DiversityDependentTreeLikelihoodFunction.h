#ifndef RB_DD_TREE_LIKELIHOOD_FUNCTION_H
#define RB_DD_TREE_LIKELIHOOD_FUNCTION_H
#include "TypedFunction.h"
#include "DiversityDependentFbdLikelihood.h"
#include "Tree.h"
namespace RevBayesCore {
template<class T> class TypedDagNode;
class DiversityDependentTreeLikelihoodFunction : public TypedFunction<double> {
public:
    DiversityDependentTreeLikelihoodFunction(const TypedDagNode<Tree>* tree,
        const TypedDagNode<double>* lambda, const TypedDagNode<double>* mu,
        const TypedDagNode<double>* alpha, const TypedDagNode<double>* capacity,
        const TypedDagNode<double>* origin, DDFBD::BirthModel rateModel,
        bool crown, unsigned condition, bool branchingTimes, DDFBD::Settings settings);
    DiversityDependentTreeLikelihoodFunction* clone() const override;
    void update() override;
protected:
    void swapParameterInternal(const DagNode*, const DagNode*) override;
private:
    const TypedDagNode<Tree>* tree;
    const TypedDagNode<double> *lambda, *mu, *alpha, *capacity, *origin;
    DDFBD::BirthModel rateModel;
    bool crown;
    unsigned condition;
    bool branchingTimes;
    DDFBD::Settings settings;
};
}
#endif
