#include "Func_speciesObservationTree.h"

#include "SpeciesObservationTreeFunction.h"
#include "RlTimeTree.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"
#include "RlString.h"
#include "RealPos.h"
#include "ModelVector.h"

using namespace RevLanguage;

/** default constructor */
Func_speciesObservationTree::Func_speciesObservationTree( void ) : TypedFunction<TimeTree>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_speciesObservationTree* Func_speciesObservationTree::clone( void ) const
{
    
    return new Func_speciesObservationTree( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::Tree>* Func_speciesObservationTree::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const TimeTree&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::SpeciesObservationTreeFunction* f = new RevBayesCore::SpeciesObservationTreeFunction( tau,
        static_cast<const ModelVector<RlString>&>(args[1].getVariable()->getRevObject()).getValue(),
        static_cast<const ModelVector<RlString>&>(args[2].getVariable()->getRevObject()).getValue(),
        static_cast<const ModelVector<RealPos>&>(args[3].getVariable()->getRevObject()).getDagNode() );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_speciesObservationTree::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "tree", TimeTree::getClassTypeSpec(), "The tree variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule("sampleNames", ModelVector<RlString>::getClassTypeSpec(), "Unique character-sample names.", ArgumentRule::BY_VALUE, ArgumentRule::CONSTANT) );
        argumentRules.push_back( new ArgumentRule("speciesNames", ModelVector<RlString>::getClassTypeSpec(), "Species endpoint name for each sample.", ArgumentRule::BY_VALUE, ArgumentRule::CONSTANT) );
        argumentRules.push_back( new ArgumentRule("ages", ModelVector<RealPos>::getClassTypeSpec(), "Actual ages of character observations.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY) );
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_speciesObservationTree::getClassType(void)
{
    
    static std::string rev_type = "Func_speciesObservationTree";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_speciesObservationTree::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}



/**
 * Get the primary Rev name for this function.
 */
std::string Func_speciesObservationTree::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnSpeciesObservationTree";
    
    return f_name;
}


const TypeSpec& Func_speciesObservationTree::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
