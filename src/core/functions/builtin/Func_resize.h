/**
 * @file
 * This file contains the declaration of Func_resize, which
 * is used to return a matrix/vector of new dimensions.
 *
 * @brief Declaration of Func_resize
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id$
 */

#ifndef Func_resize_H
#define Func_resize_H

#include "RbFunction.h"

#include <map>
#include <string>
#include <vector>

class DAGNode;
class VectorString;

const std::string Func_resize_name = "Resize function";

class Func_resize :  public RbFunction {
	
public:
	// Basic utility functions
	Func_resize*                clone(void) const;                                  //!< Clone the object
	const VectorString&         getClass(void) const;                               //!< Get class vector
    const TypeSpec&             getTypeSpec(void) const;                            //!< Get language type of the object
	
	// Regular functions
    RbPtr<const ArgumentRules>  getArgumentRules(void) const;                       //!< Get argument rules
	const TypeSpec&             getReturnType(void) const;                          //!< Get type of return value
  
protected:
	RbPtr<RbLanguageObject>     executeFunction( void);                              //!< Execute operation
  
private:
  void                          resizeVector(Container &vec, Environment *args, unsigned int numArg);
  static const TypeSpec         typeSpec;
};

#endif

#include "Ellipsis.h"
#include "RbUtil.h"
#include "TypeSpec.h"
#include "ValueRule.h"


// Definition of the static type spec member
const TypeSpec Func_resize::typeSpec(Func_resize_name);


/** Clone object */
Func_resize* Func_resize::clone( void ) const {
  
  return new Func_resize( *this );
}

void Func_resize::resizeVector(Container &vec, Environment *args, unsigned int numArg) { 
    unsigned int nrows = ( static_cast<Natural*>( (RbObject*)(*args)[numArg]->getValue() ) )->getValue();
    for (unsigned int j = 0 ; j < vec.size() ; j ++) {
        std::cout << "j: "<<j <<std::endl;
//        std::cout << "element j: "<<vec[j] <<std::endl;

        if (vec.getElement(j) != NULL && vec.getElement(j)->isTypeSpec(Vector_name )) { //if element j already a vector
            static_cast<Vector*>( (RbObject*)vec.getElement(j))->resize(nrows);
        }
        else {
            if (vec.getElement(j) != NULL) { //if element j not a vector but not NULL
                RbPtr<Vector> v2( new Vector(vec.getElement(j)->getTypeSpec()) );
                v2->push_back(vec.getElement(j));
                v2->resize(nrows);
                vec.setElement(j, RbPtr<RbObject>( (Vector*)v2));
            }
            else { //if element j NULL
                RbPtr<Vector> v2( new Vector(RbObject_name) );
                for ( size_t i = 0; i < nrows; i++ )
                    v2->push_back( RbPtr<RbObject>::getNullPtr() );
                vec.setElement(j, RbPtr<RbObject>( (Vector*)v2));
            }
        }
        if (args->size() != numArg +1) { //if we still have dimensions to populate
            resizeVector(*(static_cast<Vector*>( (RbObject*)vec.getElement(j))), args, numArg+1);
        }
        
    }
}

/** Execute function: We rely on operator overloading to provide the necessary functionality */
RbPtr<RbLanguageObject> Func_resize::executeFunction( void ) {

  unsigned long nargs = args->size();
  //The identity of the object to resize
  RbPtr<Container> val( static_cast<Container*>( (RbLanguageObject*)(*args)[0]->getValue() ) );

    if (nargs == 2) {
      //Resizing a vector
        std::cout << "nargs = 2" <<std::endl;
      unsigned int nrows = ( static_cast<Natural*>( (RbObject*)(*args)[1]->getValue() ) )->getValue();
      val->resize(nrows);
    }
    else {
        std::cout << "nargs = 3" <<std::endl;
      //Resizing a matrix of nargs-1 dimensions
       // for (unsigned int i = 0 ; i < nargs-1 ; i ++) {
            int nrows = ( static_cast<Natural*>( (RbObject*)(*args)[1]->getValue() ) )->getValue();
        
        if (!val->isType(Vector_name)) {
            val = RbPtr<Container>( static_cast<Container*>( val->convertTo(TypeSpec(Vector_name, RbPtr<TypeSpec>(new TypeSpec(RbLanguageObject_name) ) ) ) ) );
        }
        
            val->resize(nrows);
            resizeVector(*val, args, 2);
//                if (val->getElement(0)->getTypeSpec() == Vector_name) {
//                    for (unsigned int j = 0 ; j < val->size() ; j ++) {
//                        //val->setElement(j, Vector(val->getElement(j))->resize(nrows);
//                        static_cast<Vector*>( (RbObject*)val->getElement(j))->resize(nrows);
//                    }
//                }
//                else {
//                    resizeVector(val, args, 1);
// 
//                }
           // }
       // }
    }
  return RbPtr<RbLanguageObject>( val );
  
}


/** Get argument rules */
RbPtr<const ArgumentRules> Func_resize::getArgumentRules( void ) const {
  
  static RbPtr<ArgumentRules> argumentRules( new ArgumentRules() );
  static bool          rulesSet = false;
  
  if ( !rulesSet ) 
  {
    argumentRules->push_back( RbPtr<ArgumentRule>( new ValueRule( "", TypeSpec(Container_name) ) ) );
    argumentRules->push_back( RbPtr<ArgumentRule>( new ValueRule( "", Natural_name ) ) );
    argumentRules->push_back( RbPtr<ArgumentRule>( new Ellipsis( Natural_name ) ) );
    rulesSet = true;
  }
  
  return RbPtr<const ArgumentRules>( argumentRules );
}


/** Get class vector describing type of object */
const VectorString& Func_resize::getClass( void ) const {
  
  static std::string  rbName  = "Func_resize"; 
  static VectorString rbClass = VectorString( rbName ) + RbFunction::getClass();
  
  return rbClass;
}


/** Get return type */
const TypeSpec& Func_resize::getReturnType( void ) const {
	static TypeSpec rt = TypeSpec(Vector_name, RbPtr<TypeSpec>(new TypeSpec(RbLanguageObject_name) ) );
  return rt;
}


/** Get type spec */
const TypeSpec& Func_resize::getTypeSpec( void ) const {
  
  return typeSpec;
}

