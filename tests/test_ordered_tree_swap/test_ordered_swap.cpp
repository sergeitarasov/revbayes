#include "OrderedTreeSwapProposal.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "TreeChangeEventHandler.h"
#include "TreeChangeEventListener.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
using namespace RevBayesCore;

static TopologyNode* tip(const std::string& name)
{ TopologyNode* n = new TopologyNode(name); n->setAge(0.0, false); return n; }
static TopologyNode* join(TopologyNode* a, TopologyNode* b, double age)
{
    TopologyNode* n = new TopologyNode(); n->setAge(age, false);
    n->addChild(a); a->setParent(n); n->addChild(b); b->setParent(n); return n;
}
static Tree example()
{
    Tree tree; tree.setRooted(true);
    tree.setRoot(join(join(tip("A"), tip("B"), 5.0), join(tip("C"), tip("D"), 4.0), 8.0), true);
    return tree;
}
static std::string ordered(const TopologyNode& node)
{
    std::ostringstream out;
    out << node.getIndex() << ':' << node.getName() << ':' << node.getAge() << '[';
    for (const TopologyNode* child : node.getChildren()) out << ordered(*child) << ',';
    out << ']'; return out.str();
}
static size_t index(const Tree& t, const std::string& name)
{ for (const TopologyNode* n : t.getNodes()) if (n->getName() == name) return n->getIndex(); assert(false); return 0; }
class ToyDistribution : public TypedDistribution<Tree> {
public:
    explicit ToyDistribution(const Tree& t) : TypedDistribution<Tree>(new Tree(t)) {}
    ToyDistribution* clone() const { return new ToyDistribution(*this); }
    double computeLnProbability() { return 0.0; }
    void redrawValue() {}
protected:
    void swapParameterInternal(const DagNode*, const DagNode*) {}
};
class CountingListener : public TreeChangeEventListener {
public:
    unsigned count = 0;
    void fireTreeChangeEvent(const TopologyNode&, const unsigned&) { ++count; }
};
int main()
{
    GLOBAL_RNG->setSeed(7731);
    Tree t = example();
    const std::string original = ordered(t.getRoot());
    const size_t a = index(t,"A"), b = index(t,"B"), c = index(t,"C"), d = index(t,"D");
    assert(OrderedTreeSwapProposal::swapNodes(t,a,b));
    assert(t.getNode(a).getParent().getChild(0).getName() == "B");
    assert(ordered(t.getRoot()) != original);
    assert(OrderedTreeSwapProposal::swapNodes(t,a,b));
    assert(ordered(t.getRoot()) == original);
    CountingListener listener;
    t.getTreeChangeEventHandler().addListener(&listener);
    assert(OrderedTreeSwapProposal::swapNodes(t,a,d));
    assert(OrderedTreeSwapProposal::swapNodes(t,a,d));
    assert(listener.count > 0);
    t.getTreeChangeEventHandler().removeListener(&listener);
    assert(OrderedTreeSwapProposal::swapNodes(t,a,d));
    assert(&t.getNode(a).getParent() == &t.getNode(c).getParent());
    assert(OrderedTreeSwapProposal::swapNodes(t,a,d));
    assert(ordered(t.getRoot()) == original);
    assert(!OrderedTreeSwapProposal::swapNodes(t,a,t.getNode(a).getParent().getIndex()));
    // Left internal age 5 cannot move under right internal age 4.
    assert(!OrderedTreeSwapProposal::swapNodes(t,t.getNode(a).getParent().getIndex(),c));
    assert(ordered(t.getRoot()) == original);
    // Exhaustively check the self-inverse deterministic kernel and exact ordering.
    for (size_t i=0; i<t.getNumberOfNodes(); ++i) for (size_t j=0; j<t.getNumberOfNodes(); ++j)
    {
        Tree copy = t;
        const bool moved = OrderedTreeSwapProposal::swapNodes(copy,i,j);
        if (moved) assert(OrderedTreeSwapProposal::swapNodes(copy,i,j));
        assert(ordered(copy.getRoot()) == original);
        for (const TopologyNode* n : copy.getNodes())
            if (!n->isRoot()) assert(n->getBranchLength() > 0.0);
    }
    StochasticNode<Tree>* node = new StochasticNode<Tree>("tree", new ToyDistribution(t));
    OrderedTreeSwapProposal proposal(node);
    unsigned accepted = 0;
    for (unsigned rep=0; rep<100; ++rep)
    {
        proposal.prepareProposal();
        double logH = proposal.doProposal();
        if (std::isfinite(logH)) { assert(logH == 0.0); ++accepted; }
        proposal.undoProposal();
        assert(ordered(node->getValue().getRoot()) == original);
        proposal.cleanProposal();
    }
    assert(accepted > 0);
    std::cout << "Ordered tree swap reversibility tests passed.\n";
}
