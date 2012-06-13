/**
 * @file
 * This file contains the declaration of Simplex, a complex type
 * used to hold a simplex.
 *
 * @brief Declaration of Simplex
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id$
 */

#ifndef Simplex_H
#define Simplex_H

#include "MemberObject.h"
#include "RlVector.h"

#include <ostream>
#include <string>
#include <vector>


class Simplex : public MemberObject {

public:
                                Simplex(const size_t n = 2);                                //!< Simplex of length (size) n
                                Simplex(const std::vector<double>& x);                      //!< Simplex from double vector
                                Simplex(const RlVector<RealPos>& x);                        //!< Simplex from double RlVector

    double                      operator[](size_t i);                                       //!< Index op
    double                      operator[](size_t i) const;                                 //!< Const index op

    // Basic utility functions
    Simplex*                    clone(void) const;                                          //!< Clone object
    static const std::string&   getClassName(void);                                         //!< Get class name
    static const TypeSpec&      getClassTypeSpec(void);                                     //!< Get class type spec
    RbValue<void*>              getLeanValue(void) const;                                   //!< Transform the object into a basic element pointer for fast access.
    const TypeSpec&             getTypeSpec(void) const;                                    //!< Get language type of the object
    size_t                      memorySize() const;                                         //!< Get the size
    void                        printValue(std::ostream& o) const;                          //!< Print value (for user)
    void                        setLeanValue(const RbValue<void*> &val);                    //!< Set the value of the object from a basic (lean) element pointer.

    // Vector functions, including STL-like functions
    const std::vector<double>&  getValue( void ) const;
    void                        setElement(size_t index, double x);                         //!< Set the value of position
    void                        setValue(const std::vector<double>& x);                     //!< Set the value using STL vector of double
    size_t                      size(void) const;
    
private:
    void                        rescale(void);                                              //!< Rescale the simplex

    std::vector<double>         elements;
};

#endif