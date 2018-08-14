/**
 * @file DistributionBetaBinomial
 * This file contains the functions of the beta-binomial distribution.
 *
 * @brief Implementation of the beta-binomial distribution.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes core development team (JK, WD, WP)
 * @license GPL version 3
 * @version 1.0
 * @since 2011-03-17, version 1.0
 *
 * $Id$
 */


#ifndef DistributionBetaBinomial_H
#define DistributionBetaBinomial_H

namespace RevBayesCore {

    class RandomNumberGenerator;

    namespace RbStatistics {
    
        namespace BetaBinomial {
        
            double                      pdf(double y, double n, double a, double b, bool log);     /*!< Beta Binomial probability density */
            double                      lnPdf(double n, double a, double b, double y);                       /*!< Beta Binomial log_e probability density */
            double                      cdf(double n, double a, double b, double x);                                    /*!< Beta Binomial cumulative probability */
            double                      quantile(double p, double n, double a, double b);                               		    /*!< Beta Binomial(n,p) quantile */
            int                         rv(double n, double a, double b, RandomNumberGenerator& rng);             /*!< Beta Binomial random variable */
	
            double                      do_search(double y, double *z, double a, double b, double n, double pr, double incr); // ?????
        }
    }
}

#endif
