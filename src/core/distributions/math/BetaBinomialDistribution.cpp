/*
 * BetaBinomialDistribution.cpp
 *
 *  Created on: Nov 15, 2017
 *      Author: JMK
 */

#include "BetaBinomialDistribution.h"
#include "DistributionBetaBinomial.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "RbException.h"

using namespace RevBayesCore;

BetaBinomialDistribution::BetaBinomialDistribution(const TypedDagNode<long> *m, const TypedDagNode<double> *alpha, const TypedDagNode<double> *beta) : TypedDistribution<long>( new long( 0 ) ),
    n( m ),
	a( alpha ),
	b( beta )
{

    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( n );
    addParameter( a );
    addParameter( b );

    *value = RbStatistics::BetaBinomial::rv(n->getValue(), alpha->getValue(), beta->getValue(), *GLOBAL_RNG);
}


BetaBinomialDistribution::~BetaBinomialDistribution(void) {

    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



BetaBinomialDistribution* BetaBinomialDistribution::clone( void ) const {

    return new BetaBinomialDistribution( *this );
}


double BetaBinomialDistribution::computeLnProbability( void )
{

    // check that the value is inside the boundaries
    if ( *value > n->getValue() || *value < 0 )
    {
        return RbConstants::Double::neginf;
    }

    return RbStatistics::BetaBinomial::lnPdf(n->getValue(), a->getValue(), b->getValue(), *value);
}



void BetaBinomialDistribution::redrawValue( void ) {

    *value = RbStatistics::BetaBinomial::rv(n->getValue(), a->getValue(), b->getValue(), *GLOBAL_RNG);
}


/** Swap a parameter of the distribution */
void BetaBinomialDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == a)
    {
        a = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == b)
    {
        b = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == n)
    {
        n = static_cast<const TypedDagNode<long>* >( newP );
    }

}
