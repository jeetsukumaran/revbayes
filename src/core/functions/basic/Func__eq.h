/**
 * @file
 * This file contains the declaration and implementation
 * of the templated Func__eq, which is used to compare
 * two variables for equality.
 *
 * @brief Declaration and implementof Func__eq
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id$
 */

#ifndef Func__eq_H
#define Func__eq_H

#include "RbFunction.h"

#include <map>
#include <string>

class DAGNode;
class VectorString;

template <typename firstValType, typename secondValType>
class Func__eq :  public RbFunction {

    public:
        // Basic utility functions
        Func__eq*                   clone(void) const;                                          //!< Clone the object
        const VectorString&         getClass(void) const;                                       //!< Get class vector
        const TypeSpec&             getTypeSpec(void) const;                                    //!< Get language type of the object

        // Regular functions
    	RbPtr<RbLanguageObject>     execute(void);                                              //!< Execute function
        RbPtr<const ArgumentRules>  getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&             getReturnType(void) const;                                  //!< Get type of return value
    
    private:
        static const TypeSpec       typeSpec;
        static const TypeSpec       returnTypeSpec;
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
const TypeSpec Func__eq<firstValType, secondValType>::typeSpec("Func__eq", new TypeSpec(firstValType().getType() + "," + secondValType().getType()));
template <typename firstValType, typename secondValType>
const TypeSpec Func__eq<firstValType, secondValType>::returnTypeSpec(RbBoolean_name);


/** Clone object */
template <typename firstValType, typename secondValType>
Func__eq<firstValType, secondValType>* Func__eq<firstValType, secondValType>::clone( void ) const {

    return new Func__eq( *this );
}


/** Execute function: We rely on operator overloading to provide the functionality */
template <typename firstValType, typename secondValType>
RbPtr<RbLanguageObject> Func__eq<firstValType,secondValType>::execute( void ) {

    const RbPtr<firstValType>  val1( static_cast<firstValType*> ( (RbLanguageObject*)(*args)[0]->getValue() ) );
    const RbPtr<secondValType> val2( static_cast<secondValType*>( (RbLanguageObject*)(*args)[1]->getValue() ) );
    
    return RbPtr<RbLanguageObject>( new RbBoolean( *val1 == *val2 ) );
}


/** Get argument rules */
template <typename firstValType, typename secondValType>
RbPtr<const ArgumentRules> Func__eq<firstValType, secondValType>::getArgumentRules(void) const {

    static RbPtr<ArgumentRules> argumentRules( new ArgumentRules() );
    static bool          rulesSet = false;

    if ( !rulesSet ) {
        argumentRules->push_back( RbPtr<ArgumentRule>( new ValueRule( "", firstValType() .getTypeSpec() ) ) );
        argumentRules->push_back( RbPtr<ArgumentRule>( new ValueRule( "", secondValType().getTypeSpec() ) ) );
        rulesSet = true;
    }

    return RbPtr<const ArgumentRules>( argumentRules );
}


/** Get class vector describing type of object */
template <typename firstValType, typename secondValType>
const VectorString& Func__eq<firstValType, secondValType>::getClass( void ) const {

    static std::string  rbName  = "Func__eq<" + firstValType().getType() + "," + secondValType().getType() + ">"; 
    static VectorString rbClass = VectorString( rbName ) + RbFunction::getClass();
    
    return rbClass;
}


/** Get return type */
template <typename firstValType, typename secondValType>
const TypeSpec& Func__eq<firstValType, secondValType>::getReturnType( void ) const {

    return returnTypeSpec;
}


/** Get return spec */
template <typename firstValType, typename secondValType>
const TypeSpec& Func__eq<firstValType, secondValType>::getTypeSpec( void ) const {
    
    return typeSpec;
}

