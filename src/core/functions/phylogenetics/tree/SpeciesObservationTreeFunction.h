#ifndef SpeciesObservationTreeFunction_H
#define SpeciesObservationTreeFunction_H

#include "RbVector.h"
#include "Tree.h"
#include "TypedFunction.h"
#include <string>
#include <vector>

namespace RevBayesCore {
template <class valueType> class TypedDagNode;
/** Project an oriented budding-species tree onto observations at their true ages.
 * Child 0 continues the ancestral species; endpoint names identify species.
 * The initial species may have observations on the stem above the tree root.
 * Its origin constraint belongs to the upstream species-tree distribution.
 */
class SpeciesObservationTreeFunction : public TypedFunction<Tree> {
public:
    SpeciesObservationTreeFunction(const TypedDagNode<Tree>* tree,
        const std::vector<std::string>& sampleNames,
        const std::vector<std::string>& speciesNames,
        const TypedDagNode<RbVector<double> >* ages);
    SpeciesObservationTreeFunction* clone(void) const;
    void update(void);
    static Tree project(const Tree& tree, const std::vector<std::string>& sampleNames,
        const std::vector<std::string>& speciesNames, const std::vector<double>& ages);
protected:
    void swapParameterInternal(const DagNode* oldP, const DagNode* newP);
private:
    const TypedDagNode<Tree>* speciesTree;
    const TypedDagNode<RbVector<double> >* sampleAges;
    std::vector<std::string> names;
    std::vector<std::string> species;
};
}
#endif
