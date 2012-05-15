/**
 * @file
 * This file contains the implementation of Distribution, which specifies the
 * interface for Distributions in RevBayes. A Distribution is taken in
 * the statistical sense.
 *
 * @brief Declaration of Distribution
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-03-19 09:25:05 +0100 (Mon, 19 Mar 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-04-16, version 1.0
 * @interface Distribution
 *
 * $Id: Distribution.h 1353 2012-03-19 08:25:05Z hoehna $
 */

#include "Distribution.h"

#include <typeinfo>
#include <sstream>
#include <string>


/**
 * Compute the joind ln pdf.
 *
 * This is the default implementation which just computes the sum of all single lnPdf's.
 */
double Distribution::jointLnPdf( void ) const {

    double jlnpd = 0;
    
    double* singlePds = lnPdf();
    
    // compute the number of values
    const std::vector<size_t> &lengths = members[members.size()-1].lengths;
    size_t elements = 1;
    for ( size_t i = 0; i < lengths.size(); ++i) {
        elements *= lengths[i];
    }
    
    // sum up all single lnpdfs
    for ( size_t i = 0; i < elements; ++i) {
        jlnpd += singlePds[i];
    }
    
    return jlnpd;
}


/** Execute function */
double* Distribution::lnPdf( void ) const {
    
    // set all the offset and the level to 0 and delegate the function call
    std::vector<size_t> offsets;
    for (size_t i = 0; i < members.size(); ++i) {
        offsets.push_back( 0 );
    }
    
    // we allocate the memory for the lnPd's
    // compute the number of elements we'll compute
    size_t element = 1;
    size_t i = members.size()-1;
    for (size_t j = 0; j < members[i].lengths.size(); ++j) {
        element *= members[i].lengths[j];
    }
    // allocate the memory
    double* val = new double[element];
    
    // delegate the function execution
    lnPdf(val,offsets,0);
    
    // return the pointer to the lnPd's
    return val;
}


void Distribution::lnPdf(double *val, std::vector<size_t> &offset, size_t level) const {
    
    // first, we test if we have another dimension
    if ( members[members.size()-1].lengths.size() == level ) {
        // no, we execute the actual function here
        *val = lnPdfSingleValue(offset);
    }
    else {
        // next, we compute the number of elements we need to compute
        std::vector<size_t> elements;
        for (size_t i = 0; i < members.size(); ++i) {
            size_t element = 1;
            for (size_t j = level+1; j < members[i].lengths.size(); ++j) {
                element *= members[i].lengths[j];
            }
            elements.push_back( element );
        }
        
        // the number of elements for the return val
        size_t valElements = elements[elements.size()-1];
        
        // now, we iterate over the elements in this dimension
        for (size_t i = 0; i < members[members.size()-1].lengths[level]; ++i) {
            // we recompute the offset
            std::vector<size_t> nested_offsets;
            for (size_t j = 0; j < members.size(); ++j) {
                size_t nested_offset = offset[j];
                // test whether we have ele
                if ( level < members[j].lengths.size() ) {
                    nested_offset += elements[j] * i;
                }
                nested_offsets.push_back( nested_offset );
            }
            
            // and call the function for one level deeper
            lnPdf(val + (valElements*i),nested_offsets, level+1);
        }
    }
}


/** Execute function */
double* Distribution::pdf( void ) const {
    
    // set all the offset and the level to 0 and delegate the function call
    std::vector<size_t> offsets;
    for (size_t i = 0; i < members.size(); ++i) {
        offsets.push_back( 0 );
    }
    
    // we allocate the memory for the lnPd's
    // compute the number of elements we'll compute
    size_t element = 1;
    size_t i = members.size()-1;
    for (size_t j = 0; j < members[i].lengths.size(); ++j) {
        element *= members[i].lengths[j];
    }
    // allocate the memory
    double* val = new double[element];
    
    // delegate the function execution
    pdf(val,offsets,0);
    
    // return the pointer to the lnPd's
    return val;
}


void Distribution::pdf(double *val, std::vector<size_t> &offset, size_t level) const {
    
    // first, we test if we have another dimension
    if ( members[members.size()-1].lengths.size() == level ) {
        // no, we execute the actual function here
        *val = pdfSingleValue(offset);
    }
    else {
        // next, we compute the number of elements we need to compute
        std::vector<size_t> elements;
        for (size_t i = 0; i < members.size(); ++i) {
            size_t element = 1;
            for (size_t j = level+1; j < members[i].lengths.size(); ++j) {
                element *= members[i].lengths[j];
            }
            elements.push_back( element );
        }
        
        // the number of elements for the return val
        size_t valElements = elements[elements.size()-1];
        
        // now, we iterate over the elements in this dimension
        for (size_t i = 0; i < members[members.size()-1].lengths[level]; ++i) {
            // we recompute the offset
            std::vector<size_t> nested_offsets;
            for (size_t j = 0; j < members.size(); ++j) {
                size_t nested_offset = offset[j];
                // test whether we have ele
                if ( level < members[j].lengths.size() ) {
                    nested_offset += elements[j] * i;
                }
                nested_offsets.push_back( nested_offset );
            }
            
            // and call the function for one level deeper
            pdf(val + (valElements*i),nested_offsets, level+1);
        }
    }
}


/** Draw a random variable */
void Distribution::rv( void ) {
    
    // set all the offset and the level to 0 and delegate the function call
    std::vector<size_t> offsets;
    for (size_t i = 0; i < members.size(); ++i) {
        offsets.push_back( 0 );
    }
    
    // delegate the function execution
    rv(offsets,0);
    
}


void Distribution::rv( std::vector<size_t> &offset, size_t level ) {
    
    // first, we test if we have another dimension
    if ( members[members.size()-1].lengths.size() == level ) {
        // no, we execute the actual function here
        rvSingleValue(offset);
    }
    else {
        // next, we compute the number of elements we need to compute
        std::vector<size_t> elements;
        for (size_t i = 0; i < members.size(); ++i) {
            size_t element = 1;
            for (size_t j = level+1; j < members[i].lengths.size(); ++j) {
                element *= members[i].lengths[j];
            }
            elements.push_back( element );
        }
        
        // now, we iterate over the elements in this dimension
        for (size_t i = 0; i < members[members.size()-1].lengths[level]; ++i) {
            // we recompute the offset
            std::vector<size_t> nested_offsets;
            for (size_t j = 0; j < members.size(); ++j) {
                size_t nested_offset = offset[j];
                // test whether we have ele
                if ( level < members[j].lengths.size() ) {
                    nested_offset += elements[j] * i;
                }
                nested_offsets.push_back( nested_offset );
            }
            
            // and call the function for one level deeper
            rv(nested_offsets, level+1);
        }
    }
}

/** We catch here the setting of the argument variables to store our parameters. */
void Distribution::setParameters(const std::vector<RbValue<void*> > &args) {
    
    members = args;
    
    // delegate to derived class' implementation
    setInternalParameters( args );
}


/** We catch here the setting of the member variables to store our parameters. */
void Distribution::setObservedValue(const RbValue<void *> &v) {
    
    members[members.size()-1] = v;
    
    // delegate to derived class' implementation
    setInternalObservedValue( v );
    
}


std::string Distribution::toString( void ) const {
    std::string name = "";
    std::stringstream o;
    o << typeid(*this).name();
    o << " (";
    
    for (std::vector<RbValue<void*> >::const_iterator i = members.begin(); i != members.end(); ++i) {
        if ( i != members.begin() ) {
            o << ", ";
        }
        const RbValue<void*>& val = *i;
        if ( val.value == NULL ) {
            o << "NULL";
        }
        else {
            o << long( val.value );
        }
    }
    o << ")";
    
    name = o.str();
    
    return name;
}