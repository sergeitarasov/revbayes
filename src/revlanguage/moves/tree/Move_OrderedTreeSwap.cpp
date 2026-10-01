#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "OrderedTreeSwapProposal.h"
#include "RevObject.h"
#include "RealPos.h"
#include "MetropolisHastingsMove.h"
#include "Move_OrderedTreeSwap.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"

namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

Move_OrderedTreeSwap::Move_OrderedTreeSwap() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_OrderedTreeSwap* Move_OrderedTreeSwap::clone(void) const
{
    
    return new Move_OrderedTreeSwap(*this);
}


void Move_OrderedTreeSwap::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );

    RevBayesCore::Proposal *p = new RevBayesCore::OrderedTreeSwapProposal(t);
    value = new RevBayesCore::MetropolisHastingsMove(p,w);
    
}


/** Get class name of object */
const std::string& Move_OrderedTreeSwap::getClassName(void)
{
    
    static std::string rbClassName = "Move_OrderedTreeSwap";
    
    return rbClassName;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_OrderedTreeSwap::getClassTypeSpec(void)
{
    
    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rbClass;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_OrderedTreeSwap::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "OrderedTreeSwap";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_OrderedTreeSwap::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
    
        move_member_rules.push_back( new ArgumentRule( "tree", TimeTree::getClassTypeSpec(), "The tree variable this move operates on.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        
        /* Inherit weight (but not tuneTarget!) from Move and put it after the arguments created above */
        const MemberRules& inheritedRules = Move::getParameterRules();
        for (size_t i = 0; i < inheritedRules.size(); ++i)
        {
            if ( inheritedRules[i].getArgumentLabel() == "weight" )
            {
                move_member_rules.push_back( inheritedRules[i].clone() );
            }
        }
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_OrderedTreeSwap::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



/** Get type spec */
void Move_OrderedTreeSwap::printValue(std::ostream &o) const
{
    
    o << "OrderedTreeSwap(";
    if (tree != NULL)
    {
        o << tree->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a NearestNeighborInterchange variable */
void Move_OrderedTreeSwap::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}
