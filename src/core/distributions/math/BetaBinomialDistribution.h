/*
 * BetaBinomialDistribution.h
 *
 *  Created on: Nov 15, 2017
 *      Author: JMK
 */

#ifndef CORE_DISTRIBUTIONS_MATH_BETABINOMIALDISTRIBUTION_H_
#define CORE_DISTRIBUTIONS_MATH_BETABINOMIALDISTRIBUTION_H_





#endif /* CORE_DISTRIBUTIONS_MATH_BETABINOMIALDISTRIBUTION_H_ */



#ifndef BetaBinomialDistribution_H
#define BetaBinomialDistribution_H

#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {

    /**
     * @brief Beta Binomial distribution class.
     *
     * The Binomial distribution represents a family of distributions defined
     * on the natural number. The Beta Binomial distribution has 2???????? parameters:
     *   n .. the number of trials
     *   p .. the probability of success
     * Instances of this class can be associated to stochastic variables.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (JK, WD, WP)
     * @since 2013-04-12, version 1.0
     *
     */
    class BetaBinomialDistribution : public TypedDistribution<long> {

    public:
        BetaBinomialDistribution(const TypedDagNode<long> *n, const TypedDagNode<double> *p);
        virtual                                            ~BetaBinomialDistribution(void);                                             //!< Virtual destructor

        // public member functions
        BetaBinomialDistribution*                               clone(void) const;                                                      //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);

    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);        //!< Swap a parameter

    private:

        // members
        const TypedDagNode<long>*                            n;
        const TypedDagNode<double>*                         p;
    };

}

#endif
