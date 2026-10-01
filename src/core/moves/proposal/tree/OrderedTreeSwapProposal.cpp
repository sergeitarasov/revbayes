#include "OrderedTreeSwapProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "TopologyNode.h"
#include "Tree.h"
#include <algorithm>
#include <vector>

using namespace RevBayesCore;

OrderedTreeSwapProposal::OrderedTreeSwapProposal(StochasticNode<Tree>* t) :
    Proposal(), tree(t), changed(false), firstIndex(0), secondIndex(0)
{
    addNode(tree);
}
OrderedTreeSwapProposal* OrderedTreeSwapProposal::clone(void) const
{ return new OrderedTreeSwapProposal(*this); }
void OrderedTreeSwapProposal::cleanProposal(void) { changed = false; }
void OrderedTreeSwapProposal::prepareProposal(void) { changed = false; }
const std::string& OrderedTreeSwapProposal::getProposalName(void) const
{ static std::string name = "OrderedTreeSwap"; return name; }
double OrderedTreeSwapProposal::getProposalTuningParameter(void) const
{ return RbConstants::Double::nan; }
void OrderedTreeSwapProposal::printParameterSummary(std::ostream&, bool) const {}
void OrderedTreeSwapProposal::swapNodeInternal(DagNode* oldN, DagNode* newN)
{ if (oldN == tree) tree = static_cast<StochasticNode<Tree>*>(newN); }

bool OrderedTreeSwapProposal::swapNodes(Tree& tree, size_t first, size_t second)
{
    if (first == second || first >= tree.getNumberOfNodes() || second >= tree.getNumberOfNodes()) return false;
    TopologyNode* a = &tree.getNode(first);
    TopologyNode* b = &tree.getNode(second);
    if (a->isRoot() || b->isRoot()) return false;
    for (const TopologyNode* p = a; !p->isRoot(); )
    { p = &p->getParent(); if (p == b) return false; }
    for (const TopologyNode* p = b; !p->isRoot(); )
    { p = &p->getParent(); if (p == a) return false; }
    TopologyNode* pa = &a->getParent();
    TopologyNode* pb = &b->getParent();
    // Strict ordering gives valid positive-duration species branches.
    if (!(pa->getAge() > b->getAge()) || !(pb->getAge() > a->getAge()) ||
        !(pa->getAge() > a->getAge()) || !(pb->getAge() > b->getAge())) return false;
    std::vector<TopologyNode*> ca = pa->getChildren();
    std::vector<TopologyNode*> cb = pb->getChildren();
    auto ia = std::find(ca.begin(), ca.end(), a);
    auto ib = std::find(cb.begin(), cb.end(), b);
    if (ia == ca.end() || ib == cb.end()) throw RbException("OrderedTreeSwap: inconsistent parent links.");
    if (pa == pb)
    {
        const size_t sa = static_cast<size_t>(ia - ca.begin());
        const size_t sb = static_cast<size_t>(ib - cb.begin());
        std::swap(ca[sa], ca[sb]);
        std::vector<TopologyNode*> old = pa->getChildren();
        for (TopologyNode* child : old) pa->removeChild(child);
        for (TopologyNode* child : ca) pa->addChild(child);
    }
    else
    {
        *ia = b; *ib = a;
        // Rebuild exact left-to-right child arrays; addChild's position argument
        // counts from the end and ordinary remove/add would reverse orientations.
        std::vector<TopologyNode*> oldA = pa->getChildren();
        std::vector<TopologyNode*> oldB = pb->getChildren();
        for (TopologyNode* child : oldA) pa->removeChild(child);
        for (TopologyNode* child : oldB) pb->removeChild(child);
        for (TopologyNode* child : ca) pa->addChild(child);
        for (TopologyNode* child : cb) pb->addChild(child);
        a->setParent(pb);
        b->setParent(pa);
    }
    // addChild/removeChild notify topology listeners; setParent updates lengths.
    // Node indices and the root remain unchanged.
    return true;
}

double OrderedTreeSwapProposal::doProposal(void)
{
    changed = false;
    Tree& value = tree->getValue();
    std::vector<size_t> candidates;
    for (const TopologyNode* node : value.getNodes())
        if (!node->isRoot()) candidates.push_back(node->getIndex());
    if (candidates.size() < 2) return RbConstants::Double::neginf;
    size_t first = static_cast<size_t>(GLOBAL_RNG->uniform01() * candidates.size());
    size_t second = static_cast<size_t>(GLOBAL_RNG->uniform01() * (candidates.size() - 1));
    if (second >= first) ++second;
    firstIndex = candidates[first]; secondIndex = candidates[second];
    // Draw once. Resampling until valid would introduce state-dependent normalization.
    changed = swapNodes(value, firstIndex, secondIndex);
    return changed ? 0.0 : RbConstants::Double::neginf;
}
void OrderedTreeSwapProposal::undoProposal(void)
{
    if (changed)
    {
        if (!swapNodes(tree->getValue(), firstIndex, secondIndex))
            throw RbException("OrderedTreeSwap: failed to restore a valid proposal.");
        changed = false;
    }
}
