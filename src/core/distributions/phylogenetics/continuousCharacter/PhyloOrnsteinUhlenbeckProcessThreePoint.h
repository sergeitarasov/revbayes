#ifndef PhyloOrnsteinUhlenbeckProcessThreePoint_H
#define PhyloOrnsteinUhlenbeckProcessThreePoint_H

#include "AbstractPhyloContinuousCharacterProcess.h"
#include "MatrixReal.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Homogeneous distribution of character state evolution along a tree class (PhyloCTMC).
     *
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-01-23, version 1.0
     */
    class PhyloOrnsteinUhlenbeckProcessThreePoint : public AbstractPhyloContinuousCharacterProcess {
        
    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        PhyloOrnsteinUhlenbeckProcessThreePoint(const TypedDagNode<Tree> *t, size_t ns );
        virtual                                                            ~PhyloOrnsteinUhlenbeckProcessThreePoint(void);                                                              //!< Virtual destructor
        
        // public member functions
        // pure virtual
        virtual PhyloOrnsteinUhlenbeckProcessThreePoint*                    clone(void) const;                                                                      //!< Create an independent clone
        
        // non-virtual
        double                                                              computeLnProbability(void);
        void                                                                setAlpha(const TypedDagNode< double >* a);
        void                                                                setRootState(const TypedDagNode< double >* s);
        void                                                                setRootState(const TypedDagNode< RbVector< double > >* s);
        void                                                                setSigma(const TypedDagNode< double >* s);
        void                                                                setTheta(const TypedDagNode< double >* t);
        
        
    protected:
        // virtual methods that may be overwritten, but then the derived class should call this methods
        void                                                                recursivelyFlagNodeDirty(const TopologyNode& n);
        void                                                                resetValue( void );
        void                                                                simulateRecursively(const TopologyNode& node, std::vector< ContinuousTaxonData > &t);
        std::vector<double>                                                 simulateRootCharacters(size_t n);
        double                                                              sumRootLikelihood(void);
        virtual void                                                        touchSpecialization(DagNode *toucher, bool touchAll);
        
        // Parameter management functions.
        virtual void                                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                         //!< Swap a parameter
        
    private:
        double                                                              computeRootState(size_t siteIdx) const;
        void                                                                recursiveComputeRootToTipDistance( std::vector<double> &m, double v, const TopologyNode &n, size_t ni );
        std::set<size_t>                                                    recursiveComputeDistanceMatrix( MatrixReal &m, const TopologyNode &node, size_t node_index );
        
        const TypedDagNode< double >*                                       alpha;
        const TypedDagNode< double >*                                       root_state;
        const TypedDagNode< double >*                                       sigma;
        const TypedDagNode< double >*                                       theta;
        const TypedDagNode< RbVector< double > >*                           site_specific_root_state;
        
        size_t                                                              num_tips;
        std::vector<std::vector<double> >                                   obs;
        
    };
    
}


#endif
