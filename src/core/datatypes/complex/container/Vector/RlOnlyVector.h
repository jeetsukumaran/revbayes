/**
 * @file
 * This file contains the declaration of Vector, a container type
 * that acts as a base class for all constant Vectors.
 *
 * This class is a wrapper for the stl-class Vector and we provide additional RevBayes functionality
 * (e.g. getTypeSpec and getClassTypeSpec for argument checking). Furthermore, this class can only be used
 * for the RevLanguage environment. It does not provide an interface to the RevBayes core!
 *
 *
 * @brief Declaration of Vector
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-06-01 16:05:37 +0200 (Fri, 01 Jun 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-12-04, version 1.0
 *
 * $Id: RlOnlyVector.h 1605 2012-06-01 14:05:37Z hoehna $
 */

#ifndef RlOnlyVector_H
#define RlOnlyVector_H

#include "Container.h"
#include "MethodTable.h"
#include "SimpleMemberFunction.h"

#include <iostream>
#include <vector>


template <typename rlType>
class RlOnlyVector : public Container {
    
public:
    typedef typename std::vector<rlType*>::iterator iterator;
    typedef typename std::vector<rlType*>::const_iterator const_iterator;
       
    RlOnlyVector(void);                                                                                             //!< Default constructor with type RbLanguageObject
    RlOnlyVector(const std::vector<rlType*> &v);                                                                    //!< Default constructor with type RbLanguageObject
    RlOnlyVector(const rlType &v);                                                                                  //!< Default constructor with type RbLanguageObject
    RlOnlyVector(size_t n);                                                                                         //!< Default constructor with type RbLanguageObject
    RlOnlyVector(size_t n, const rlType &v);                                                                        //!< Default constructor with type RbLanguageObject
    RlOnlyVector(const RlOnlyVector& v);                                                                            //!< Copy Constructor
    
    virtual                                        ~RlOnlyVector(void);                                             //!< Virtual destructor 
    
    // Basic utility functions 
    RlOnlyVector*                                   clone(void) const;                                              //!< Clone object
    RbObject*                                       convertTo(const TypeSpec& type) const;                          //!< Convert to type
    static const std::string&                       getClassName(void);                                             //!< Get class name
    static const TypeSpec&                          getClassTypeSpec(void);                                         //!< Get class type spec
    const TypeSpec&                                 getTypeSpec(void) const;                                        //!< Get language type of the object
    virtual bool                                    isConvertibleTo(const TypeSpec& type) const;                    //!< Is this object convertible to the asked one?
    virtual void                                    printValue(std::ostream& o) const;                              //!< Print value for user
    
    rlType&                                         operator[](size_t index);                                       //!< subscript operator
    const rlType&                                   operator[](size_t index) const;                                 //!< subscript operator (const)
    RlOnlyVector&                                   operator=(const RlOnlyVector& x);                               //!< Assignment operator
    RlOnlyVector&                                   operator+=(const rlType& x);                                    //!< Concatenate
    RlOnlyVector&                                   operator+=(const RlOnlyVector& x);                              //!< Concatenate
    const RlOnlyVector                              operator+(const rlType& x) const;                               //!< Concatenate
    const RlOnlyVector                              operator+(const RlOnlyVector& x) const;                         //!< Concatenate
    bool                                            operator==(const RlOnlyVector& x) const;                        //!< Equality
    bool                                            operator!=(const RlOnlyVector& x) const;                        //!< Inequality
        
    
    // Member object function
    RbPtr<RbLanguageObject>                         executeSimpleMethod(const std::string& name, const std::vector<const RbObject*>& args);         //!< Override to map member methods to internal functions
    const MemberRules&                              getMemberRules(void) const;                                     //!< Get member rules
    const MethodTable&                              getMethods(void) const;                                         //!< Get methods
    
    // Vector functions
    iterator                                        begin(void);                                                    //!< Iterator to the beginning of the Vector
    const_iterator                                  begin(void) const;                                              //!< Const-iterator to the beginning of the Vector
    void                                            clear(void);                                                    //!< Clear
    iterator                                        end(void);                                                      //!< Iterator to the end of the Vector
    const_iterator                                  end(void) const;                                                //!< Const-iterator to the end of the Vector
    int                                             findIndex(const RbObject& x) const;                             //!< Find the index if the element being equal to that one
    RbPtr<RbObject>                                 getElement(size_t index);                                       //!< Get element (non-const to return non-const element)
//    const std::vector<rlType*>&                     getValue(void) const;                                           //!< Get the stl Vector of elements
    void                                            pop_back(void);                                                 //!< Drop element at back
    void                                            pop_front(void);                                                //!< Drop element from front
    void                                            push_back(const rlType &x);                                     //!< Append element to end
    void                                            push_front(const rlType &x);                                    //!< Append element to end
    void                                            resize(size_t n);                                               //!< Resize to new AbstractVector of length n
    void                                            setElement(const size_t index, RbObject *elem);                 //!< Set element with type conversion
    size_t                                          size(void) const;                                               //!< get the number of elements in the AbstractVector
    
protected:
    
    // We store internally pointers to our objects. This is necessary because elements can be also of the derived type and we need to be able to make proper copies of the Vector and all its elements
    std::vector<rlType*>                            elements;
    
private:
    
    MemberRules                                     memberRules;
    MethodTable                                     methods;
    
};



#include "ArgumentRule.h"
#include "Complex.h"
#include "Ellipsis.h"
#include "MethodTable.h"
#include "Monitor.h"
#include "Move.h"
#include "Probability.h"
#include "RbBoolean.h"
#include "RbException.h"
#include "RbNullObject.h"
#include "RbString.h"
#include "RbUtil.h"
#include "SimpleMemberFunction.h"
#include "TypeSpec.h"

#include <algorithm>

/** Vector type of elements */
template <typename rlType>
RlOnlyVector<rlType>::RlOnlyVector( void ) : Container( rlType::getClassTypeSpec() ) {
    
}


/** Constructor with dimension (n) and NULL pointers to every object */
template <typename rlType>
RlOnlyVector<rlType>::RlOnlyVector(size_t n) : Container( rlType::getClassTypeSpec() )  {
    
    for (size_t i = 0; i < n; i++) {
        this->push_back( NULL );
    }
}


/** Constructor with dimension (n) and copys of x for every object */
template <typename rlType>
RlOnlyVector<rlType>::RlOnlyVector(size_t n, const rlType &x) : Container( rlType::getClassTypeSpec() )  {
    
    for (size_t i = 0; i < n; i++) {
        this->push_back( x.clone() );
    }
}


/** Constructor with dimension (n) and copys of x for every object */
template <typename rlType>
RlOnlyVector<rlType>::RlOnlyVector(const std::vector<rlType*> &x) : Container( rlType::getClassTypeSpec() )  {
    
    elements = x;
}



/** Copy Constructor */
template <typename rlType>
RlOnlyVector<rlType>::RlOnlyVector(const RlOnlyVector<rlType> &v) : Container(v), memberRules( v.memberRules ), methods( v.methods ) {
    
    typename std::vector<rlType*>::const_iterator it;
    // copy all the elements by deep copy
    for ( it = v.elements.begin(); it != v.elements.end(); it++) {
        elements.push_back( (*it)->clone() );
    }
    
}


/** Destructor. Free the memory of the elements. */
template <typename rlType>
RlOnlyVector<rlType>::~RlOnlyVector(void) {
    
    // just call clear which will free the memory of the elements
    clear();
}

/** Assignment operator; make sure we get independent elements */
template <typename rlType>
RlOnlyVector<rlType>& RlOnlyVector<rlType>::operator=( const RlOnlyVector<rlType>& x ) {
    
    if ( this != &x ) {
        
        // First assign using parent assignment operator. This will test to make sure the containers
        // are of the same type, and throw an error if they are not. By calling it before we destroy
        // our own elements, we can make sure that an assignment error leaves us intact, which it should
        Container::operator=( x );
        
        // just call clear which will free the memory of the objects
        clear();
        
        typename std::vector<typename rlType::valueType>::const_iterator i;
        for ( i = x.elements.begin(); i != x.elements.end(); i++ ) {
            elements.push_back( (*i)->clone() );
        }
        
        memberRules     = x.memberRules;
        methods         = x.methods;
    }
    
    return ( *this );
}


/* Subscript operator */
template <typename rlType>
rlType& RlOnlyVector<rlType>::operator[]( size_t index ) {
    
    return *elements[index];
}


/* Subscript operator */
template <typename rlType>
const rlType& RlOnlyVector<rlType>::operator[]( size_t index ) const {
    
    return *elements[index];
}


/** Concatenation with operator+ (valueType) */
template <typename rlType>
RlOnlyVector<rlType>& RlOnlyVector<rlType>::operator+=( const rlType& x ) {
    
    push_back( x.clone() );
    
    return *this;
}


/** Concatenation with operator+ (RlOnlyVector) */
template <typename rlType>
RlOnlyVector<rlType>& RlOnlyVector<rlType>::operator+=( const RlOnlyVector<rlType>& x ) {
    
    for ( size_t i = 0; i < x.elements.size(); i++ )
        push_back( x[i].clone() );
    
    return *this;
}


/** Equals comparison */
template <typename rlType>
bool RlOnlyVector<rlType>::operator==(const RlOnlyVector<rlType>& x) const {
    
    if (size() != x.size())
        return false;
    
    for (size_t i=0; i<elements.size(); i++) {
        if (elements[i] != x[i])
            return false;
    }
    
    return Container::operator==( x );
}


/** Not equals comparison */
template <typename rlType>
bool RlOnlyVector<rlType>::operator!=(const RlOnlyVector<rlType>& x) const {
    
    return !operator==(x);
}


/** Concatenation with operator+ (valueType) */
template <typename rlType>
const RlOnlyVector<rlType> RlOnlyVector<rlType>::operator+( const rlType& x ) const {
    
    RlOnlyVector tempVec = *this;
    
    tempVec.push_back( x.getValue() );
    
    return tempVec;
}


/** Concatenation with operator+ (RlOnlyVector) */
template <typename rlType>
const RlOnlyVector<rlType> RlOnlyVector<rlType>::operator+( const RlOnlyVector<rlType>& x ) const {
    
    RlOnlyVector<rlType> tempVec = *this;
    
    for ( size_t i = 0; i < x.elements.size(); i++ )
        tempVec.push_back( x[i] );
    
    return tempVec;
}


/** Get iterator to the beginning of the Vector. */
template <typename rlType>
typename std::vector<rlType*>::iterator RlOnlyVector<rlType>::begin( void ) {
    return elements.begin();
}


/** Get iterator to the beginning of the Vector. */
template <typename rlType>
typename std::vector<rlType*>::const_iterator RlOnlyVector<rlType>::begin( void ) const {
    return elements.begin();
}


/** Convertible to: default implementation */
template <typename rlType>
RbObject* RlOnlyVector<rlType>::convertTo(const TypeSpec &type) const {
    
    // test whether we want to convert to another Vector
    if ( type.getBaseType() == getClassName() ) {
        
        // work through all the possible base element types
        
        // RbLanguageObject
//        if ( type.getElementType() == Integer::getClassTypeSpec() ) {
//            RlOnlyVector<Integer>* convObject = new RlOnlyVector<Integer>();
//            // insert copies of all objects. clone if they are of the right type, otherwise convert them
//            typename std::vector<rlType*>::const_iterator i;
//            for ( i = begin(); i != end(); i++) {
//                rlType orgElement = rlType( *i );
//                // test whether this element is already of the right type
//                if ( orgElement.isTypeSpec(type.getElementType()) ) {
//                    convObject->push_back( dynamic_cast<Integer &>( orgElement ) );
//                }
//                else {
//                    convObject->push_back( *static_cast<Integer*>( orgElement.convertTo(type.getElementType()) ) );
//                }
//            }
//            
//            return convObject;
//        }
        
        
    }
    
    return Container::convertTo(type);
}



/** Clear contents of value container */
template <typename rlType>
void RlOnlyVector<rlType>::clear( void ) {
    
    typename std::vector<rlType*>::iterator i;
    for ( i = elements.begin(); i != elements.end(); i++) {
        RbObject* theElement = *i;
        delete theElement;
    }
    
    elements.clear();
}


template <typename rlType>
RlOnlyVector<rlType>* RlOnlyVector<rlType>::clone() const {
    return new RlOnlyVector<rlType>( *this );
}


/** Get iterator to the end of the Vector. */
template <typename rlType>
typename std::vector<rlType*>::iterator RlOnlyVector<rlType>::end( void ) {
    return elements.end();
}


/** Get iterator to the end of the Vector. */
template <typename rlType>
typename std::vector<rlType*>::const_iterator RlOnlyVector<rlType>::end( void ) const {
    return elements.end();
}


/** Execute member function. */
template <typename rlType>
RbPtr<RbLanguageObject> RlOnlyVector<rlType>::executeSimpleMethod(std::string const &name, const std::vector<const RbObject *> &args) {
    
//    if ( name == "sort" ) {
//        sort();
//        
//        return NULL;
//    }
//    else if ( name == "unique" ) {
//        unique();
//        
//        return NULL;
//    }
    
    return Container::executeSimpleMethod(name, args);
}



/**
 * Find the index of the given element.
 * We rely on overloaded operator== in the element classes to check for matches.
 * 
 * \param x the element we are looking for. 
 * \return The index or -1 if we didn't find it.
 */
template <typename rlType>
int RlOnlyVector<rlType>::findIndex(const RbObject& x) const {
    
    // get the iterator to the first element
    typename std::vector<typename rlType::valueType>::const_iterator i;
    
    // initialize the index
    int index = 0;
    for ( i = elements.begin(); i != elements.end(); i++, index++) {
        const RbObject& element = *(*i);
        
        // test if the object matches
        // note that we rely on the implemented operator==
        if ( element == x ) {
            return index;
        }
    }
    
    return -1;
}


/* Get class name of object */
template <typename rlType>
const std::string& RlOnlyVector<rlType>::getClassName(void) { 
    
    static std::string rbClassName = "Vector";
    
	return rbClassName; 
}

/* Get class type spec describing type of object */
template <typename rlType>
const TypeSpec& RlOnlyVector<rlType>::getClassTypeSpec(void) { 
    
    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Container::getClassTypeSpec() ), new TypeSpec( rlType::getClassTypeSpec() ) );
    
	return rbClass; 
}


/* Get element */
template <typename rlType>
RbPtr<RbObject> RlOnlyVector<rlType>::getElement(size_t index) {
    
    return elements[index];
}


/** 
 * Get the member rules
 * We expect that a Vector is created by "Vector(x,...)". 
 * All variables are for simplicity just single elements. For more sophisticated constructors (e.g. from a vector of elements)
 * constructor functions should be used.
 */
template <typename rlType>
const MemberRules& RlOnlyVector<rlType>::getMemberRules(void) const {
    
    static MemberRules memberRules;
    static bool rulesSet = false;
    
    if ( ! rulesSet ) {
        // set the member rules
        memberRules.push_back( new ArgumentRule( "x", true, rlType::getClassTypeSpec() ) );
        memberRules.push_back( new Ellipsis( rlType::getClassTypeSpec() ) );
        
        rulesSet = true;
    }
    
    return memberRules;
}


/** Get the methods for this vector class */
/* Get method specifications */
template <typename rlType>
const MethodTable& RlOnlyVector<rlType>::getMethods(void) const {
    
    static MethodTable methods;
    static bool methodsSet = false;
    
    if (!methodsSet) {
        
        // add method for call "x[]" as a function
        ArgumentRules* squareBracketArgRules = new ArgumentRules();
        squareBracketArgRules->push_back( new ArgumentRule( "index" , true, Natural::getClassTypeSpec() ) );
        methods.addFunction("[]",  new SimpleMemberFunction( rlType::getClassTypeSpec(), squareBracketArgRules) );
        
//        // add method for call "x.sort()" as a function
//        ArgumentRules* sortArgRules = new ArgumentRules();
//        methods.addFunction("sort",  new SimpleMemberFunction( RbVoid_name, sortArgRules) );
//        
//        // add method for call "x.unique()" as a function
//        ArgumentRules* uniqueArgRules = new ArgumentRules();
//        methods.addFunction("unique",  new SimpleMemberFunction( RbVoid_name, uniqueArgRules) );
        
        // necessary call for proper inheritance
        methods.setParentTable( &Container::getMethods() );
        
        methodsSet = true;
    }
    
    
    return methods;
}


/** Get the type spec of this class. We return a member variable because instances might have different element types. */
template <typename rlType>
const TypeSpec& RlOnlyVector<rlType>::getTypeSpec(void) const {
    
    static TypeSpec typeSpec = getClassTypeSpec();
    return typeSpec;
}


/* Is convertible to: default implementation */
template <typename rlType>
bool RlOnlyVector<rlType>::isConvertibleTo(const TypeSpec &type) const {
    
    // test whether we want to convert to another Vector
    if ( type.getBaseType() == getClassName() ) {
        
        // work through all the possible base element types
        typename std::vector<rlType*>::const_iterator i;
        for ( i = begin(); i != end(); i++) {
            rlType* orgElement = *i;
            // test whether this element is already of the right type
            if ( !orgElement->isTypeSpec(type.getElementType()) && !orgElement->isConvertibleTo(type.getElementType()) ) {
                return false;
            }
        }
        
        return true;
    }
    
    return false;
}


/** Print value for user */
template <typename rlType>
void RlOnlyVector<rlType>::printValue( std::ostream& o ) const {
    
    o << "[ ";
    typename std::vector<rlType*>::const_iterator i;
    for ( i = elements.begin(); i != elements.end(); i++ ) {
        if ( i != elements.begin() )
            o << ", ";
        rlType* tmp = *i;
        tmp->printValue(o);
    }
    o <<  " ]";
    
}


/** Pop element off of front of vector, updating length in process */
template <typename rlType>
void RlOnlyVector<rlType>::pop_front(void) {
    
    elements.erase(elements.begin());
}


/** Pop element off of back of vector, updating length in process */
template <typename rlType>
void RlOnlyVector<rlType>::pop_back(void) {
    
    elements.pop_back();
}


/** Push an int onto the back of the vector */
template <typename rlType>
void RlOnlyVector<rlType>::push_back( const rlType &x ) {
    
    elements.push_back( x.clone() );
    
}


/** Push an int onto the front of the vector */
template <typename rlType>
void RlOnlyVector<rlType>::push_front( const rlType &x ) {
    
    elements.insert( elements.begin(), x.clone() );
}


/** Resize vector */
template <typename rlType>
void RlOnlyVector<rlType>::resize( size_t n ) {
    
    if ( n < elements.size() )
        throw RbException( "Invalid attempt to shrink vector" );
    
    for ( size_t i = elements.size(); i < n; i++ )
        elements.push_back( NULL );
}


/* Set element */
template <typename rlType>
void RlOnlyVector<rlType>::setElement(const size_t index, RbObject *elem) {
    if (index >= elements.size()) {
        throw RbException("Cannot set element in Vector outside the current range.");
    }
    
    // remove first the old element at the index
    elements.erase(elements.begin()+index);
    
    throw RbException("Missing implementation of RlOnlyVector::setElement()");
    //    elements.insert(elements.begin()+index, *elem);
}


/** Get the size of the vector */
template <typename rlType>
size_t RlOnlyVector<rlType>::size( void ) const {
    
    return elements.size();
    
}



#endif