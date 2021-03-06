#ifndef NegativeBinomialDistribution_H
#define NegativeBinomialDistribution_H

#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {
    
    /**
     * @brief Binomial distribution class.
     *
     * The Binomial distribution represents a family of distributions defined
     * on the natural number. The Binomial distribution has 2 parameters:
     *   n .. the number of trials
     *   p .. the probability of success
     * Instances of this class can be associated to stochastic variables.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2013-04-12, version 1.0
     *
     */
    class NegativeBinomialDistribution : public TypedDistribution<long> {
        
    public:
        NegativeBinomialDistribution(const TypedDagNode<long> *r, const TypedDagNode<double> *p);
        virtual                                            ~NegativeBinomialDistribution(void);                                             //!< Virtual destructor
        
        // public member functions
        NegativeBinomialDistribution*                       clone(void) const;                                                      //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<long>*                           r;
        const TypedDagNode<double>*                         p;
    };
    
}

#endif
