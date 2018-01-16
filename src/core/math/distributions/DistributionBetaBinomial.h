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
        
            double                      pdf(double n, double p, double x);                                    /*!< Binomial(n,p) probability density */
            double                      pdf(double y, double n, double pp, double a, double b, bool log);     /*!< Beta Binomial(n,pp) probability density */
            double                      lnPdf(double n, double pp, double a, double b, double y);                       /*!< Beta Binomial(n,pp) log_e probability density */
            double                      cdf(double n, double p, double x);                                    /*!< Beta Binomial(n,p) cumulative probability */
            double                      quantile(double n, double p);                               		    /*!< Beta Binomial(n,p) quantile */
            int                         rv(double n, double pp, double a, double b, RandomNumberGenerator& rng);             /*!< Beta Binomial(n,pp) random variable */
	
            double                      do_search(double y, double *z, double p, double n, double pr, double incr); // ?????
        }
    }
}

#endif
