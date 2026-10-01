#include "DiversityDependentTreeLikelihoodFunction.h"
#include "TypedDagNode.h"
#include "TopologyNode.h"
#include "RbException.h"
#include <algorithm>
#include <cmath>
#include <functional>
using namespace RevBayesCore;
DiversityDependentTreeLikelihoodFunction::DiversityDependentTreeLikelihoodFunction(
    const TypedDagNode<Tree>* t, const TypedDagNode<double>* l,
    const TypedDagNode<double>* m, const TypedDagNode<double>* a,
    const TypedDagNode<double>* c, const TypedDagNode<double>* o,
    DDFBD::BirthModel rate, bool cr, unsigned co, bool bt, DDFBD::Settings se)
    : TypedFunction<double>(new double(0)), tree(t), lambda(l), mu(m), alpha(a),
      capacity(c), origin(o), rateModel(rate), crown(cr), condition(co),
      branchingTimes(bt), settings(se)
{
    addParameter(tree); addParameter(lambda); addParameter(mu);
    addParameter(alpha); addParameter(capacity); addParameter(origin);
    update();
}
DiversityDependentTreeLikelihoodFunction* DiversityDependentTreeLikelihoodFunction::clone() const
{ return new DiversityDependentTreeLikelihoodFunction(*this); }
void DiversityDependentTreeLikelihoodFunction::swapParameterInternal(const DagNode* oldP, const DagNode* newP)
{
    if (oldP == tree) tree = static_cast<const TypedDagNode<Tree>*>(newP);
    if (oldP == lambda) lambda = static_cast<const TypedDagNode<double>*>(newP);
    if (oldP == mu) mu = static_cast<const TypedDagNode<double>*>(newP);
    if (oldP == alpha) alpha = static_cast<const TypedDagNode<double>*>(newP);
    if (oldP == capacity) capacity = static_cast<const TypedDagNode<double>*>(newP);
    if (oldP == origin) origin = static_cast<const TypedDagNode<double>*>(newP);
}
void DiversityDependentTreeLikelihoodFunction::update()
{
    const Tree& t = tree->getValue();
    std::vector<double> ages;
    for (const TopologyNode* node : t.getNodes()) {
        if (!std::isfinite(node->getAge()) || node->getAge() < 0 ||
            (!node->isRoot() && node->getAge() > node->getParent().getAge()))
            throw RbException(RbException::MATH_ERROR, "DD tree has invalid ages");
        if (node->isTip()) {
            if (node->getAge() != 0) throw RbException("DD tree likelihood requires extant tips at age zero");
        } else {
            if (node->getNumberOfChildren() != 2) throw RbException("DD tree likelihood requires a binary tree");
            ages.push_back(node->getAge());
        }
    }
    std::sort(ages.begin(),ages.end(),std::greater<double>());
    DDFBD::Parameters p{lambda->getValue(),alpha->getValue(),mu->getValue(),0,1,rateModel,capacity->getValue()};
    try {
        *value = DDFBD::extantTreeLogLikelihood(crown ? t.getRoot().getAge() : origin->getValue(),
            ages,crown,condition,branchingTimes,p,settings);
    } catch (const std::exception& e) {
        throw RbException(RbException::MATH_ERROR,e.what());
    }
}
