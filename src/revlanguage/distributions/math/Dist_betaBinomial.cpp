#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "BetaBinomialDistribution.h"
#include "ContinuousStochasticNode.h"
#include "Dist_betaBinomial.h"
#include "Probability.h"
#include "RealPos.h"
#include "RbException.h"
#include "DistributionBetaBinomial.h"
#include "RbHelpReference.h"

using namespace RevLanguage;

Dist_betaBinomial::Dist_betaBinomial(void) : TypedDistribution<Natural>()
{

}


Dist_betaBinomial::~Dist_betaBinomial(void)
{

}



Dist_betaBinomial* Dist_betaBinomial::clone( void ) const
{

    return new Dist_betaBinomial(*this);
}


RevBayesCore::BetaBinomialDistribution* Dist_betaBinomial::createDistribution( void ) const
{

    // get the parameters
    RevBayesCore::TypedDagNode<long>*    vn = static_cast<const Natural     &>( n->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* va = static_cast<const RealPos &>( a->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* vb = static_cast<const RealPos &>( b->getRevObject() ).getDagNode();
    RevBayesCore::BetaBinomialDistribution* d  = new RevBayesCore::BetaBinomialDistribution( vn, va, vb );
    return d;
}



/* Get Rev type of object */
const std::string& Dist_betaBinomial::getClassType(void)
{

    static std::string rev_type = "Dist_betaBinomial";
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_betaBinomial::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
	return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_betaBinomial::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "betaBinomial";

    return d_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_betaBinomial::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );

    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::vector<std::string> Dist_betaBinomial::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::vector<std::string> descriptions;
    descriptions.push_back( "Beta Binomial probability distribution of x successes in n trials." );

    return descriptions;
}


/**
 * Get the more detailed description of the function
 */
std::vector<std::string> Dist_betaBinomial::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::vector<std::string> details;

    std::string details_1 = "";
    details_1 += "The binomial probability distribution defines the number of success in n trials,";
    details_1 += "where each trial has the same success probability p. The probability is given by";
    details_1 += "(n choose x) p^(x) * (1-p)^(n-p)";

    details.push_back( details_1 );

    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_betaBinomial::getHelpExample(void) const //FIX THIS
{
    // create an example as a single string variable.
    std::string example = "";

    example += "p ~ dnBeta(1.0,1.0)\n";
    example += "x ~ dnBinomial(n=10,p)\n";
    example += "x.clamp(8)\n";
    example += "moves[1] = mvSlide(p, delta=0.1, weight=1.0)\n";
    example += "monitors[1] = screenmonitor(printgen=1000, separator = \"\t\", x)\n";
    example += "mymodel = model(p)\n";
    example += "mymcmc = mcmc(mymodel, monitors, moves)\n";
    example += "mymcmc.burnin(generations=20000,tuningInterval=100)\n";
    example += "mymcmc.run(generations=200000)\n";

    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_betaBinomial::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;


    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_betaBinomial::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnBernoulli" );


    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_betaBinomial::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Beta Binomial Distribution";

    return title;
}


/** Return member rules (no members) */
const MemberRules& Dist_betaBinomial::getParameterRules(void) const
{

    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( !rules_set )
    {

        dist_member_rules.push_back( new ArgumentRule( "n", Natural::getClassTypeSpec()    , "Number of trials.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "a", RealPos::getClassTypeSpec()    , "alpha parameter", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "b", RealPos::getClassTypeSpec()    , "beta parameter", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }

    return dist_member_rules;
}


const TypeSpec& Dist_betaBinomial::getTypeSpec( void ) const
{

    static TypeSpec ts = getClassTypeSpec();
    return ts;
}


/** Print value for user */
void Dist_betaBinomial::printValue(std::ostream& o) const
{

    o << "BetaBinomial(p=";
    if ( a != NULL )
    {
        o << b->getName();
    }
    else
    {
        o << "?";
    }
    o << ", n=";
    if ( n != NULL )
        o << n->getName();
    else
        o << "?";
    o << ")";
}


/** Set a member variable */
void Dist_betaBinomial::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "n" )
    {
        n = var;
    }
    else if ( name == "a" )
    {
        a = var;
    }
    else if ( name == "b" )
    {
        b = var;
    }
    else
    {
        TypedDistribution<Natural>::setConstParameter(name, var);
    }
}


