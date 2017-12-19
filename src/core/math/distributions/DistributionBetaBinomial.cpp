/**
 * @file DistributionBetaBinomial
 * This file contains the functions of the beta-binomial distribution.
 *
 * @brief Implementation of the beta binomial distribution.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes core development team (JK, WD, WP) 15.11.2017
 * @license GPL version 3
 * @version 1.0
 * @since 2011-03-17, version 1.0
 *
 * $Id$
 */

#include <cmath>

#include "DistributionBeta.h"
#include "DistributionBinomial.h"
#include "DistributionNormal.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathFunctions.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbMathHelper.h"
#include "RbMathLogic.h"
#include "DistributionBetaBinomial.h"
#include <cmath>

using namespace RevBayesCore;

/*!
 * This function calculates the probability density 
 * for a beta-binomially-distributed random variable.
 * The beta-binomial distribution is the binomial distribution
 * in which the probability of successes at each trial is random,
 * and follows a beta distribution.
 *
 * \brief Beta Binomial probability density.
 * \param n is the number of trials. 
 * \param p is the success probability. 
 * \param x is the number of successes. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::BetaBinomial::cdf(double n, double p, double x)
{
	throw RbException("The Beta Binomial cdf is not yet implemented in RB.");
}

/*!
 * This function draws a random variable
 * from a beta-binomial distribution.
 *
 * From R:
 *
 *     (1) pdf() has both p and q arguments, when one may be represented
 *         more accurately than the other (in particular, in df()).
 *     (2) pdf() does NOT check that inputs x and n are integers. This
 *         should be done in the calling function, where necessary.
 *         -- but is not the case at all when called e.g., from df() or dbeta() !
 *     (3) Also does not check for 0 <= p <= 1 and 0 <= q <= 1 or NaN's.
 *         Do this in the calling function.
 *
 *  REFERENCE
 *
 *	Kachitvichyanukul, V. and Schmeiser, B. W. (1988).
 *	Binomial random variate generation.
 *	Communications of the ACM 31, 216-222.
 *	(Algorithm BTPEC).
 *
 *
 * \brief Beta-Binomial probability density.
 * \param n is the number of trials.
 * \param pp is the success probability.
 * \param a is the number of successes. ????
 * \param b is the ????
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */


int RbStatistics::BetaBinomial::rv(double n, double pp, double a, double b, RevBayesCore::RandomNumberGenerator &rng)
{
	int y;

	double p = RbStatistics::Beta::rv(a,b,rng);
	y = RbStatistics::Binomial::rv(n,pp,rng);
	return y;
}

/*!
 * This function calculates the probability density 
 * for a beta-binomially-distributed random variable.
 *
 * \brief Beta-Binomial probability density.
 * \param n is the number of trials. 
 * \param p is the success probability. 
 * \param x is the number of successes. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double RbStatistics::BetaBinomial::lnPdf(double n, double pp, double a, double b) {

    double aa = 1.0 - pp;
    return lnPdf(n, pp, aa, b); //is this right?
}

/*!
 * This function calculates the probability density 
 * for a beta-binomially-distributed random variable.
 *
 * \brief Beta-Binomial probability density.
 * \param n is the number of trials. 
 * \param pp is the success probability.
 * \param y is the number of successes.
 * \param a is the alpha parameter for the beta distribution
 * \param b is the beta parameter for the beta distribution
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */

//double RbStatistics::BetaBinomial::pdf(double y, double n, double pp, double a, double b, bool asLog)
//{

//}


double RbStatistics::BetaBinomial::pdf(double y, double n, double pp, double a, double b, bool asLog)
{

    double constant;
    if(pp == 0)
    		return((y == 0) ? (asLog ? 0.0 : 1.0 ) : (asLog ? RbConstants::Double::neginf : 0.0));

    //TODO: conditionals
    if (a==0)
    		return((y == 0) ? (asLog ? 0.0 : 1.0) : (asLog ? RbConstants::Double::neginf : 0.0) );
    		//return constant;

    if (b==0)
    		return((y == n) ? (asLog ? 0.0 : 1.0) : (asLog ? RbConstants::Double::neginf : 0.0) );
    		//return constant;

    constant = RevBayesCore::RbMath::lnChoose(n, y);
    if(asLog == false){
    		double prUnnorm = constant * RbStatistics::Beta::pdf(a-y, b+n-y, y);
    		double prNormed = prUnnorm / RbStatistics::Beta::pdf(a, b, y);
    		return prNormed;
    }
    else{
    		RbMath::lnChoose = RevBayesCore::RbMathFunctions::ln(constant); //error: no member named 'ln' in RbMathFunctions
    		double prUnnorm = RbMath::lnChoose + RbStatistics::Beta::lnPdf(a, b, y);
    		double prNormed = prUnnorm - RbStatistics::Beta::lnPdf(a, b, y);
    }



}


/*!
 * This function calculates the probability density 
 * for a beta-binomially-distributed random variable.
 *
 * From R:
 *
 *     (1) pdf() has both p and q arguments, when one may be represented
 *         more accurately than the other (in particular, in df()).
 *     (2) pdf() does NOT check that inputs x and n are integers. This
 *         should be done in the calling function, where necessary.
 *         -- but is not the case at all when called e.g., from df() or dbeta() !
 *     (3) Also does not check for 0 <= p <= 1 and 0 <= q <= 1 or NaN's.
 *         Do this in the calling function.
 *
 * \brief Beta Binomial probability density.
 * \param n is the number of trials. 
 * \param p is the success probability. 
 * \param x is the number of successes. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */

//double RbStatistics::BetaBinomial::quantile(double quantile_prob, double n, double p)
double quantile(double quantile_prob, double n, double p)
{
	throw RbException("There is no simple formula for this, and it is not yet implemented in RB.");
}
