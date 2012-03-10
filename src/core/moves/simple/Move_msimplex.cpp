/**
 * @file
 * This file contains the implementation of Move_msimplex,
 * which changes a simplex.
 *
 * @brief Implementation of Move_msimplex
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-09-08, version 1.0
 *
 * $Id$
 */

#include "Distribution.h"
#include "DistributionDirichlet.h"
#include "Move_msimplex.h"
#include "Natural.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "RbUtil.h"
#include "RbStatisticsHelper.h"
#include "RealPos.h"
#include "Simplex.h"
#include "StochasticNode.h"
#include "ValueRule.h"
#include "Workspace.h"

#include <cmath>


/** Constructor for parser */
Move_msimplex::Move_msimplex( void ) : MoveSimple( getMemberRules() ), alpha( NULL ), numCategories( NULL ) {
}


/** Clone object */
Move_msimplex* Move_msimplex::clone( void ) const {

    return new Move_msimplex( *this );
}


/** Get class name of object */
const std::string& Move_msimplex::getClassName(void) { 
    
    static std::string rbClassName = "Simplex move";
    
	return rbClassName; 
}

/** Get class type spec describing type of object */
const TypeSpec& Move_msimplex::getClassTypeSpec(void) { 
    
    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( MoveSimple::getClassTypeSpec() ) );
    
	return rbClass; 
}

/** Get type spec */
const TypeSpec& Move_msimplex::getTypeSpec( void ) const {
    
    static TypeSpec typeSpec = getClassTypeSpec();
    
    return typeSpec;
}


/** Return member rules */
const MemberRules& Move_msimplex::getMemberRules( void ) const {

    static MemberRules memberRules = MemberRules();
    static bool        rulesSet = false;

    if (!rulesSet) 
    {
        memberRules.push_back( new ValueRule( "variable", Simplex::getClassTypeSpec() ) );

        /* Inherit weight from MoveSimple, put it after variable */
        const MemberRules& inheritedRules = MoveSimple::getMemberRules();
        memberRules.insert( memberRules.end(), inheritedRules.begin(), inheritedRules.end() ); 

        memberRules.push_back( new ValueRule( "tuning"  , RealPos::getClassTypeSpec() ) );
        memberRules.push_back( new ValueRule( "num_cats", Natural::getClassTypeSpec() ) );

        rulesSet = true;
		}

    return memberRules;
}



/** Return the random variable type appropriate for the move */
const TypeSpec Move_msimplex::getVariableType( void ) const {

    return Simplex::getClassTypeSpec();
}


/** Perform the move */
double Move_msimplex::perform( void ) {

    // Get random number generator    
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    // Get relevant values
    StochasticNode*        nodePtr = static_cast<StochasticNode*>( node->getDagNode() );
    double                 alpha0  = static_cast<const RealPos&>( alpha->getValue()   ).getValue();
    int                    k       = static_cast<const Natural&>( numCategories->getValue() ).getValue();

	std::vector<double> curVal = static_cast<const Simplex&>( nodePtr->getValue() ).getValue();
	std::vector<double> newVal = curVal;
    int                 n      = int( curVal.size() );

	/* We update the simplex values by proposing new values from a Dirichlet centered
	   on the current values. The i-th parameter of the Dirichlet is the i-th value
	   in the simplex multiplied by a parameter (alpha0, AKA tuning) that controls the
	   variance of the Dirichlet. We implement two cases of this general move. In one
	   case, all of the elements of the simplex are targeted for update (n == k). In the
	   other, more complicated, case a subset of the elements of the simplex are updated
	   (k < n). Here, we construct a smaller simplex with k+1 elements. The first k of the
	   elements are the values from the full simplex that were targeted for update. The last
	   element of the smaller simplex accumulates the probabilities of all of the simplex
	   values in the full simplex that were not targeted for update. We then update the
	   small simplex by centering a Dirichlet on the small simplex. The values for those elements
	   in the full simplex that were not targeted for update are all changed proportionally.
	   This means that we need to calculate the Jacobian for the Hastings ratio in this case. */
	double lnProposalRatio = 0.0;
	if ( k > n )
		{
		// we can't update more elements than there are elements in the simplex
		throw RbException( "Attempting to update too many simplex variables" );
		}
	else if ( k < 1 )
		{
		// we also can't update 0 or a negative number of elements
		throw RbException( "Attempting to update too few simplex variables" );
		}
	else if ( k < n )
		{
		// we update a subset of the elements in the full simplex
		// pick k values at random, producing a map from the index in the full vector (curVal) to
		// the index in the reduced vector (x, alpha, and z)
		std::vector<int> indicesToUpdate;
		std::vector<int> tmpV;
		for (int i=0; i<n; i++)
			tmpV.push_back(i);
		RbStatistics::Helper::randomlySelectFromVectorWithoutReplacement<int>(tmpV, indicesToUpdate, k, *rng);
		std::map<int,int> mapper;
		for (size_t i=0; i<indicesToUpdate.size(); i++)
			mapper.insert( std::make_pair(indicesToUpdate[i], i) );
			
		// set up the vectors
		std::vector<double> x(indicesToUpdate.size()+1, 0.0);
		std::vector<double> alphaForward(indicesToUpdate.size()+1, 0.0);
		std::vector<double> alphaReverse(indicesToUpdate.size()+1, 0.0);
		std::vector<double> z(indicesToUpdate.size()+1, 0.0);
		for (int i=0; i<n; i++)
			{
			std::map<int,int>::iterator it = mapper.find(i);
			if (it != mapper.end())
				x[it->second] += curVal[it->first];
			else 
				x[x.size()-1] += curVal[i];
			}
		for (size_t i=0; i<x.size(); i++)
			alphaForward[i] = x[i] * alpha0;
			
		// draw a new value for the reduced vector
		z = RbStatistics::Dirichlet::rv( alphaForward, *rng );
		
		// fill in the Dirichlet parameters for the reverse probability calculations
		for (size_t i=0; i<z.size(); i++)
			alphaReverse[i] = z[i] * alpha0;
		
		// fill in the full vector
		double factor = z[z.size()-1] / x[x.size()-1];
		for (int i=0; i<n; i++)
			{
			std::map<int,int>::iterator it = mapper.find(i);
			if (it != mapper.end())
				newVal[i] = z[it->second];
			else 
				newVal[i] = curVal[i] * factor;
			}
			
		// Hastings ratio
		lnProposalRatio  = RbStatistics::Dirichlet::lnPdf(alphaReverse, x) - RbStatistics::Dirichlet::lnPdf(alphaForward, z); // Hastings Ratio
		lnProposalRatio += (n - k) * log(factor); // Jacobian
		}
	else 
		{
		// we update all of the elements in the vector
		// first, we get the parameters of the Dirichlet for the forward move
		std::vector<double> alphaForward(curVal.size());
		for (size_t i=0; i<curVal.size(); i++)
			alphaForward[i] = curVal[i] * alpha0;
			
		// then, we propose new values
		newVal = RbStatistics::Dirichlet::rv( alphaForward, *rng );
		
		// and calculate the Dirichlet parameters for the (imagined) reverse move
		std::vector<double> alphaReverse(newVal.size());
        for (size_t i=0; i<curVal.size(); i++) {
			alphaReverse[i] = newVal[i] * alpha0;
            // we need to check for 0 values
            if (alphaReverse[i] < 1E-10) {
                // very low proposal probability which will hopefully result into a rejected proposal
                return -1E10;
            }
        }
		
		// finally, we calculate the log of the Hastings ratio
		lnProposalRatio = RbStatistics::Dirichlet::lnPdf(alphaReverse, curVal) - RbStatistics::Dirichlet::lnPdf(alphaForward, newVal);
		}
		
    nodePtr->setValue( new Simplex( newVal ) );
	
    return lnProposalRatio;
}



/** We catch here the setting of the member variables to store our parameters. */
void Move_msimplex::setMemberVariable(std::string const &name, Variable* var) {
    
    // test whether we want to set the variable 
    if ( name == "tuning" ) {
        alpha = var;
    } 
    else if ( name == "num_cats" ) {
        numCategories = var;
    }
    else {
        MoveSimple::setMemberVariable(name, var);
    }
}

