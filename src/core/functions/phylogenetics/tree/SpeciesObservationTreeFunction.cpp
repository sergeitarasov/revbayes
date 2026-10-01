#include "SpeciesObservationTreeFunction.h"
#include "RbException.h"
#include "TopologyNode.h"
#include "TreeChangeEventHandler.h"
#include "TypedDagNode.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <set>

using namespace RevBayesCore;

SpeciesObservationTreeFunction::SpeciesObservationTreeFunction(const TypedDagNode<Tree>* t,
    const std::vector<std::string>& n, const std::vector<std::string>& s,
    const TypedDagNode<RbVector<double> >* a) : TypedFunction<Tree>(new Tree()),
    speciesTree(t), sampleAges(a), names(n), species(s)
{
    addParameter(speciesTree);
    addParameter(sampleAges);
    update();
}

SpeciesObservationTreeFunction* SpeciesObservationTreeFunction::clone(void) const
{
    return new SpeciesObservationTreeFunction(*this);
}

void SpeciesObservationTreeFunction::update(void)
{
    // Dynamic domain failures propagate as MATH_ERROR for automatic MH rejection.
    Tree projected = project(speciesTree->getValue(), names, species, sampleAges->getValue());
    *value = projected;
    // Tree assignment preserves listeners but does not emit node changes.
    // Signal all nodes because every projected topology/branch may have changed.
    for (TopologyNode* node : value->getNodes())
        value->getTreeChangeEventHandler().fire(*node);
}

void SpeciesObservationTreeFunction::swapParameterInternal(const DagNode* oldP, const DagNode* newP)
{
    if (oldP == speciesTree) speciesTree = static_cast<const TypedDagNode<Tree>*>(newP);
    if (oldP == sampleAges) sampleAges = static_cast<const TypedDagNode<RbVector<double> >*>(newP);
}

Tree SpeciesObservationTreeFunction::project(const Tree& tree,
    const std::vector<std::string>& names, const std::vector<std::string>& species,
    const std::vector<double>& ages)
{
    if (names.empty() || names.size() != species.size() || names.size() != ages.size())
        throw RbException("fnSpeciesObservationTree requires nonempty equally sized sampleNames, speciesNames and ages.");
    std::map<std::string, const TopologyNode*> endpoints;
    for (const TopologyNode* node : tree.getNodes())
    {
        if (!std::isfinite(node->getAge()) || node->getAge() < 0.0)
            throw RbException(RbException::MATH_ERROR, "fnSpeciesObservationTree: species-tree ages must be finite and nonnegative.");
        if (!node->isRoot() && node->getAge() > node->getParent().getAge())
            throw RbException(RbException::MATH_ERROR, "fnSpeciesObservationTree: species tree has a negative branch length.");
        if (!node->isTip() && node->getNumberOfChildren() != 2)
            throw RbException("fnSpeciesObservationTree requires a bifurcating oriented species tree.");
        if (node->isTip() && (node->getName().empty() || !endpoints.emplace(node->getName(), node).second))
            throw RbException("fnSpeciesObservationTree: species endpoint names must be nonempty and unique.");
    }
    std::map<const TopologyNode*, std::vector<size_t> > edgeSamples;
    std::map<std::string, size_t> sampleIndices;
    for (size_t i = 0; i < names.size(); ++i)
    {
        if (names[i].empty() || !sampleIndices.emplace(names[i], i).second)
            throw RbException("fnSpeciesObservationTree: sample names must be nonempty and unique.");
        if (!std::isfinite(ages[i]) || ages[i] < 0.0)
            throw RbException(RbException::MATH_ERROR, "fnSpeciesObservationTree: observation ages must be finite and nonnegative.");
        const auto found = endpoints.find(species[i]);
        if (found == endpoints.end())
            throw RbException("fnSpeciesObservationTree: unknown species '" + species[i] + "'.");
        const TopologyNode* node = found->second;
        if (ages[i] < node->getAge())
            throw RbException(RbException::MATH_ERROR, "fnSpeciesObservationTree: observation '" + names[i] + "' is younger than its species endpoint.");
        // Ascend only ancestral-continuation edges. Reaching child 1 is birth.
        while (!node->isRoot() && ages[i] > node->getParent().getAge())
        {
            const TopologyNode* parent = &node->getParent();
            if (&parent->getChild(0) != node)
                throw RbException(RbException::MATH_ERROR, "fnSpeciesObservationTree: observation '" + names[i] + "' predates its species birth.");
            node = parent;
        }
        edgeSamples[node].push_back(i);
    }
    for (auto& entry : edgeSamples)
        std::stable_sort(entry.second.begin(), entry.second.end(), [&](size_t a, size_t b) { return ages[a] < ages[b]; });

    // Every recursive result owns its whole subtree until transferred to a parent.
    typedef std::unique_ptr<TopologyNode> NodePtr;
    std::function<NodePtr(const TopologyNode&)> build = [&](const TopologyNode& node) -> NodePtr {
        NodePtr below;
        if (!node.isTip())
        {
            NodePtr left = build(node.getChild(0));
            NodePtr right = build(node.getChild(1));
            if (left && right)
            {
                below.reset(new TopologyNode());
                below->setAge(node.getAge(), false);
                below->addChild(left.get()); left->setParent(below.get()); left.release();
                below->addChild(right.get()); right->setParent(below.get()); right.release();
            }
            else below = left ? std::move(left) : std::move(right);
        }
        const auto samples = edgeSamples.find(&node);
        if (samples != edgeSamples.end()) for (size_t i : samples->second)
        {
            NodePtr tip(new TopologyNode());
            Taxon taxon(names[i]);
            taxon.setSpeciesName(species[i]);
            taxon.setAge(ages[i]);
            taxon.setExtinct(ages[i] > 0.0);
            tip->setTaxon(taxon);
            tip->setAge(ages[i], false);
            if (below)
            {
                NodePtr parent(new TopologyNode());
                parent->setAge(ages[i], false);
                parent->addChild(below.get()); below->setParent(parent.get()); below.release();
                parent->addChild(tip.get()); tip->setParent(parent.get());
                tip->setSampledAncestor(true); tip.release();
                below = std::move(parent);
            }
            else below = std::move(tip);
        }
        return below;
    };
    NodePtr root = build(tree.getRoot());
    Tree result;
    result.setRooted(true);
    result.setRoot(root.release(), true);
    // A stable data-tip order is essential when a topology changes during MCMC.
    size_t internalIndex = names.size();
    for (TopologyNode* node : result.getNodes())
        node->setIndex(node->isTip() ? sampleIndices.at(node->getName()) : internalIndex++);
    result.orderNodesByIndex();
    return result;
}
