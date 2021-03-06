#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_setwd.h"
#include "RbException.h"
#include "RbSettings.h"
#include "RlString.h"
#include "RlUtils.h"
#include "TypeSpec.h"

#include <fstream>

using namespace RevLanguage;

/** Default constructor */
Func_setwd::Func_setwd( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_setwd* Func_setwd::clone( void ) const
{
    
    return new Func_setwd( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_setwd::execute( void )
{
    
    const std::string &wd = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();
    
    std::cerr << "Received directory:\t\t" << wd << std::endl;
    
    RbSettings& s = RbSettings::userSettings();
    s.setWorkingDirectory( wd );
    
    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_setwd::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "wd", RlString::getClassTypeSpec(), "The new working directory.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_setwd::getClassType(void)
{
    
    static std::string rev_type = "Func_setwd";
    
	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_setwd::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_setwd::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "setwd";
    
    return f_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_setwd::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::vector<std::string> Func_setwd::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::vector<std::string> descriptions;
    descriptions.push_back( "Set the current working directory which RevBayes uses." );
    
    return descriptions;
}


/**
 * Get the more detailed description of the function
 */
std::vector<std::string> Func_setwd::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::vector<std::string> details;
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Func_setwd::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    example += "# get the current working directory\n";
    example += "getwd()\n";
    example += "\n";
    example += "# let us set a new working directory\n";
    example += "setwd(\"~/Desktop\")\n";
    example += "\n";
    example += "# check the working directory again\n";
    example += "getwd()\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Func_setwd::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_setwd::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "getwd" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Func_setwd::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Set and print the working directory";
    
    return title;
}



/** Get type spec */
const TypeSpec& Func_setwd::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_setwd::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlUtils::Void;
    
    return return_typeSpec;
}

