#ifndef Func_WattersonTheta_H
#define Func_WattersonTheta_H

#include "RealPos.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of Watterson's theta function.
     *
     * The RevLanguage wrapper of Watternson's theta function connects
     * the variables/parameters of the function and creates the internal WattersonThetaFunction object.
     * Please read the WattersonThetaFunction.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-04-30, version 1.0
     *
     */
    class Func_WattersonTheta : public TypedFunction<RealPos> {
        
    public:
        Func_WattersonTheta( void );
        
        // Basic utility functions
        Func_WattersonTheta*                        clone(void) const;                                                              //!< Clone the object
        static const std::string&                   getClassType(void);                                                             //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                         //!< Get class type spec
        const TypeSpec&                             getTypeSpec(void) const;                                                        //!< Get the type spec of the instance
        
        // Function functions you have to override
        RevBayesCore::TypedFunction< double >*      createFunction(void) const;                                                     //!< Create internal function object
        const ArgumentRules&                        getArgumentRules(void) const;                                                   //!< Get argument rules
        
    };
    
}

#endif
