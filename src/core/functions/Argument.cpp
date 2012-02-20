/**
 * @file
 * This file contains the implementation of Argument, which is
 * used to hold a potentially labeled argument passed to a
 * function.
 *
 * @brief Implementation of Argument
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-11-20, version 1.0
 *
 * $Id$
 */

#include <sstream>

#include "Argument.h"
#include "DAGNode.h"
#include "RbUtil.h"
#include "VectorString.h"


// Definition of the static type spec member
const TypeSpec Argument::typeSpec(Argument_name);

/** Construct from argument label and DAG node */
Argument::Argument(const RbVariablePtr& v) : RbInternal() {
    
    label   = "";
    var     = v;
}


/** Construct from argument label and DAG node */
Argument::Argument(const std::string& argLabel, const RbVariablePtr& v) : RbInternal() {

    label   = argLabel;
    var     = v;
}

/** Copy Constructor. We keep the same pointer to the variable stored inside this argument. */
Argument::Argument(const Argument &x) : RbInternal(x) {
    
    label   = x.label;
    if (x.var != NULL)
        var     = x.var;
}


/** Destructor */
Argument::~Argument() {
    
}


Argument& Argument::operator=(const Argument &x) {
    
    if ( &x != this ) {
        
        if (var != NULL) {
            delete var;
        }
        // Copy the new variable
        if (x.var == NULL) {
            var = NULL;
        }
        else {
            var     = x.var;
        }
        label   = x.label;
    }
    
    return (*this);
}


/** Get class vector describing type of object */
const VectorString& Argument::getClass(void) const { 

    static VectorString rbClass = VectorString(Argument_name) + RbInternal::getClass();
	return rbClass; 
}

const std::string& Argument::getLabel(void) const {
    return label;
}


/** Get the type spec of this class. We return a static class variable because all instances will be exactly from this type. */
const TypeSpec& Argument::getTypeSpec(void) const {
    return typeSpec;
}


const Variable& Argument::getVariable(void) const {
    return *var;
}


Variable& Argument::getVariable(void) {
    return *var;
}


const RbVariablePtr& Argument::getVariablePtr(void) const {
    return var;
}


/** Complete info about object */
void Argument::printValue(std::ostream &o) const {

    o << label << " = ";
    var->printValue(o);

}


/** Set the variable of the argument */
void Argument::setVariable(const RbVariablePtr& newVar) {
    
    var = newVar;
}

