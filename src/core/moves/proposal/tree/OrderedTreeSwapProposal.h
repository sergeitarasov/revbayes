#ifndef OrderedTreeSwapProposal_H
#define OrderedTreeSwapProposal_H
#include "Proposal.h"
#include <cstddef>
namespace RevBayesCore {
class Tree;
class DagNode;
template <class T> class StochasticNode;
/** Uniform nonroot-node swap preserving ordered child slots, including siblings. */
class OrderedTreeSwapProposal : public Proposal {
public:
    explicit OrderedTreeSwapProposal(StochasticNode<Tree>* tree);
    OrderedTreeSwapProposal* clone(void) const;
    void cleanProposal(void);
    double doProposal(void);
    const std::string& getProposalName(void) const;
    double getProposalTuningParameter(void) const;
    void prepareProposal(void);
    void printParameterSummary(std::ostream&, bool) const;
    void undoProposal(void);
    // Deterministic, self-inverse kernel, also exposed for exact reversibility tests.
    // Returns false without mutation for an invalid pair.
    static bool swapNodes(Tree& tree, size_t first, size_t second);
protected:
    void swapNodeInternal(DagNode* oldN, DagNode* newN);
private:
    StochasticNode<Tree>* tree;
    bool changed;
    size_t firstIndex;
    size_t secondIndex;
};
}
#endif
