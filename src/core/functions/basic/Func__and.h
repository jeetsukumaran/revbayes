/**
 * @file
 * This file contains the declaration and implementation
 * of the templated Func__and, which is used to perform
 * the logical and operation.
 *
 * @brief Declaration and implementation of Func__and
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id$
 */

#ifndef Func__and_H
#define Func__and_H

#include "RbFunction.h"
#include "RbBoolean.h"

#include <string>

class DAGNode;
class VectorString;

template <typename firstValType, typename secondValType>
class Func__and :  public RbFunction {

    public:
    
        // Basic utility functions
        Func__and*                  clone(void) const;                                          //!< Clone the object
        const VectorString&         getClass(void) const;                                       //!< Get class vector
        const TypeSpec&             getTypeSpec(void) const;                                    //!< Get language type of the object

        // Regular functions
        const ArgumentRules&        getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&             getReturnType(void) const;                                  //!< Get type of return value

    protected:
        const RbLanguageObject&     executeFunction(void);                                      //!< Execute function
        void                        setArgumentVariable(const std::string& name, const RbVariablePtr& var);

    private:
        static const TypeSpec       typeSpec;
        static const TypeSpec       returnTypeSpec;
    
        // Arguments
        RbVariablePtr               first;
        RbVariablePtr               second;
    
        // function return value
        RbBoolean                   retValue;
};

#endif


#include "RbBoolean.h"
#include "DAGNode.h"
#include "Integer.h"
#include "MatrixReal.h"
#include "RbException.h"
#include "RbUtil.h"
#include "Real.h"
#include "TypeSpec.h"
#include "ValueRule.h"
#include "VectorString.h"


// Definition of the static type spec member
template <typename firstValType, typename secondValType>
const TypeSpec Func__and<firstValType, secondValType>::typeSpec("Func__and", new TypeSpec(firstValType().getType() + "," + secondValType().getType() ) );
template <typename firstValType, typename secondValType>
const TypeSpec Func__and<firstValType, secondValType>::returnTypeSpec(RbBoolean_name);


/** Clone object */
template <typename firstValType, typename secondValType>
Func__and<firstValType, secondValType>* Func__and<firstValType, secondValType>::clone( void ) const {

    return new Func__and( *this );
}


/** Execute function: We rely on operator overloading to provide the functionality */
template <typename firstValType, typename secondValType>
const RbLanguageObject& Func__and<firstValType,secondValType>::executeFunction( void ) {

    const firstValType&  val1 = static_cast<firstValType&> ( first->getValue()  );
    const secondValType& val2 = static_cast<secondValType&>( second->getValue() );
    retValue.setValue( val1 && val2 );
    return retValue;
}


/** Get argument rules */
template <typename firstValType, typename secondValType>
const ArgumentRules& Func__and<firstValType, secondValType>::getArgumentRules( void ) const {

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rulesSet = false;

    if ( !rulesSet ) 
        {
        argumentRules.push_back( new ValueRule( "first", firstValType() .getTypeSpec() ) );
        argumentRules.push_back( new ValueRule( "second", secondValType().getTypeSpec() ) );
        rulesSet = true;
        }

    return argumentRules;
}


/** Get class vector describing type of object */
template <typename firstValType, typename secondValType>
const VectorString& Func__and<firstValType, secondValType>::getClass( void ) const {

    static std::string  rbName  = "Func__and<" + firstValType().getType() + "," + secondValType().getType() + ">"; 
    static VectorString rbClass = VectorString( rbName ) + RbFunction::getClass();
    
    return rbClass;
}


/** Get return type */
template <typename firstValType, typename secondValType>
const TypeSpec& Func__and<firstValType, secondValType>::getReturnType( void ) const {

    return returnTypeSpec;
}


/** Get return spec */
template <typename firstValType, typename secondValType>
const TypeSpec& Func__and<firstValType, secondValType>::getTypeSpec( void ) const {
    
    return typeSpec;
}


/** We catch here the setting of the argument variables to store our parameters. */
template <typename firstValType, typename secondValType>
void Func__and<firstValType, secondValType>::setArgumentVariable(std::string const &name, const RbVariablePtr& var) {
    
    if ( name == "first" ) {
        first = var;
    }
    else if ( name == "second" ) {
        second = var;
    }
    else {
        RbFunction::setArgumentVariable(name, var);
    }
}

