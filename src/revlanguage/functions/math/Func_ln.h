#ifndef Func_ln_H
#define Func_ln_H

#include "Real.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the natural logarithm function.
     *
     * The RevLanguage wrapper of the sum function connects
     * the variables/parameters of the function and creates the internal LnFunction object.
     * Please read the LnFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-07-27, version 1.0
     *
     */
    class Func_ln :  public TypedFunction<Real> {
        
    public:
        Func_ln( void );
        
        // Basic utility functions
        Func_ln*                                        clone(void) const;                                          //!< Clone the object
        static const std::string&                       getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                     getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction<double>*            createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

#endif
