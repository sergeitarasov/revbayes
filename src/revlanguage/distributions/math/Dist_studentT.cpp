//
//  Dist_studentT.cpp
//  revbayes_work
//
//  Created by Dismukes, Wade T [EEOBS] on 11/4/16.
//  Copyright © 2016 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#include "Dist_studentT.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "StudentTDistribution.h"
#include "ContinuousStochasticNode.h"
#include "Natural.h"
#include "Probability.h"
#include "RealPos.h"

using namespace RevLanguage;

Dist_studentT::Dist_studentT(void) : ContinuousDistribution()
{
    
}


Dist_studentT::~Dist_studentT(void)
{
    
}



Dist_studentT* Dist_studentT::clone( void ) const
{
    
    return new Dist_studentT(*this);
}


RevBayesCore::StudentTDistribution* Dist_studentT::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<long>*    df = static_cast<const Natural     &>( degrees->getRevObject() ).getDagNode();
    RevBayesCore::StudentTDistribution* d  = new RevBayesCore::StudentTDistribution( df );
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_studentT::getClassType(void)
{
    
    static std::string rev_type = "Dist_studentT";
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_studentT::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    return rev_type_spec;
}



std::vector<std::string> Dist_studentT::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "studentT" );
    a_names.push_back( "t" );
    a_names.push_back( "gossetT" );
    
    return a_names;
}

/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_studentT::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "StudentT";
    
    return d_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */



/**
 * Get the (brief) description for this function
 */



/**
 * Get the more detailed description of the function
 */



/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */



/**
 * Get some references/citations for this function
 *
 */



/**
 * Get the names of similar and suggested other functions
 */



/**
 * Get the title of this help entry
 */



/** Return member rules (no members) */
const MemberRules& Dist_studentT::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        dist_member_rules.push_back( new ArgumentRule( "df", Natural::getClassTypeSpec(), "The degrees of freedom.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_studentT::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    return ts;
}


/** Print value for user */
void Dist_studentT::printValue(std::ostream& o) const
{
    
    o << "StudentT(df=";
    if ( degrees != NULL )
    {
        o << degrees->getName();
    }
    else
    {
        o << "?";
    }
    
    o << ")";
}


/** Set a member variable */
void Dist_studentT::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "df" )
    {
        degrees = var;
    }
    else
    {
        TypedDistribution<Real>::setConstParameter(name, var);
    }
}
