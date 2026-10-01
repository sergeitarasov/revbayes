// Standalone core test: link with RevBayes core and SpeciesObservationTreeFunction.
#include "SpeciesObservationTreeFunction.h"
#include "TopologyNode.h"
#include "RbException.h"
#include <cassert>
#include <iostream>
using namespace RevBayesCore;

static Tree twoSpecies()
{
    TopologyNode* a = new TopologyNode("A"); a->setAge(0.0, false);
    TopologyNode* b = new TopologyNode("B"); b->setAge(1.0, false);
    TopologyNode* root = new TopologyNode(); root->setAge(5.0, false);
    root->addChild(a); a->setParent(root);
    root->addChild(b); b->setParent(root);
    Tree tree; tree.setRooted(true); tree.setRoot(root, true); return tree;
}
static void expectFailure(const Tree& tree, const std::vector<std::string>& n,
                          const std::vector<std::string>& s, const std::vector<double>& a)
{
    bool failed = false;
    try { SpeciesObservationTreeFunction::project(tree, n, s, a); }
    catch (const RbException&) { failed = true; }
    assert(failed);
}
int main()
{
    Tree species = twoSpecies();
    Tree projected = SpeciesObservationTreeFunction::project(species,
        {"A_old", "A_mid", "A_now", "B_fossil"}, {"A", "A", "A", "B"}, {7.0, 3.0, 0.0, 2.0});
    assert(projected.getNumberOfTips() == 4);
    assert(projected.getNumberOfNodes() == 7);
    assert(projected.getRoot().getAge() == 7.0);
    const std::vector<TopologyNode*>& nodes = projected.getNodes();
    assert(nodes[0]->getName() == "A_old" && nodes[0]->getAge() == 7.0);
    assert(nodes[0]->isSampledAncestorTip());
    assert(nodes[1]->getName() == "A_mid" && nodes[1]->getAge() == 3.0);
    assert(nodes[1]->isSampledAncestorTip());
    assert(nodes[2]->getAge() == 0.0 && !nodes[2]->isSampledAncestorTip());
    assert(nodes[3]->getAge() == 2.0 && !nodes[3]->isSampledAncestorTip());
    for (const TopologyNode* node : nodes) if (!node->isRoot())
        assert(node->getBranchLength() >= 0.0);
    // Pruning unsampled endpoints suppresses their bifurcation, but not observations.
    Tree serial = SpeciesObservationTreeFunction::project(species,
        {"a", "b"}, {"A", "A"}, {4.0, 2.0});
    assert(serial.getRoot().getAge() == 4.0 && serial.getNumberOfNodes() == 3);
    assert(serial.getNodes()[0]->isSampledAncestorTip());
    // A fossil-only species is observed at its fossil age, not its extinction time.
    Tree singleton = SpeciesObservationTreeFunction::project(species, {"b"}, {"B"}, {3.0});
    assert(singleton.getNumberOfNodes() == 1 && singleton.getRoot().getAge() == 3.0);
    expectFailure(species, {"b"}, {"B"}, {6.0}); // before budding birth
    expectFailure(species, {"b"}, {"B"}, {0.0}); // after extinction
    expectFailure(species, {"a"}, {"unknown"}, {2.0});
    expectFailure(species, {"a", "a"}, {"A", "B"}, {2.0, 2.0});
    // Input species tree remains untouched.
    assert(species.getNumberOfTips() == 2 && species.getRoot().getAge() == 5.0);
    std::cout << "Species observation tree projection tests passed.\n";
}
