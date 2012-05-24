/**
 * @file
 * This file contains the implementation of CharacterStateDiscrete, which is
 * the abstract base class for discrete character data types in RevBayes.
 *
 * @brief Implementation of CharacterStateDiscrete
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id$
 */

#include "CharacterStateDiscrete.h"
#include "ConstArgumentRule.h"
#include "MemberObject.h"
#include "Natural.h"
#include "RbException.h"
#include "RbString.h"
#include "RbUtil.h"
#include "RbVector.h"



/** Constructor taking the number of states */
CharacterStateDiscrete::CharacterStateDiscrete(size_t n) : Character(), numStates(n) {

    isGapState = false;
    value.resize(numStates);
}


/** Get class name of object */
const std::string& CharacterStateDiscrete::getClassName(void) { 
    
    static std::string rbClassName = "Discrete character";
    
	return rbClassName; 
}

/** Get class type spec describing type of object */
const TypeSpec& CharacterStateDiscrete::getClassTypeSpec(void) { 
    
    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Character::getClassTypeSpec() ) );
    
	return rbClass; 
}


const MemberRules& CharacterStateDiscrete::getMemberRules(void) const {
    
    static MemberRules memberRules = MemberRules();
    static bool        rulesSet = false;
    
    if (!rulesSet) {
        
        memberRules.push_back( new ConstArgumentRule( "state"  , RbString::getClassTypeSpec() ) );
       
        rulesSet = true;
    }
    
    return memberRules;
}


/** Return the number of observed states for this character (any number greater than one indicates ambiguity) */
size_t CharacterStateDiscrete::getNumberObservedStates(void) const {

    size_t nOn = 0;
    for (size_t i=0; i<value.size(); i++)
        {
        if (value[i] == true)
            nOn++;
        }
    return nOn;
}


/** Get the unsigned representation for this character */
unsigned CharacterStateDiscrete::getUnsignedValue(void) const {

    if ( sizeof(unsigned) * 8 < value.size() )
        throw ( RbException("Too many states to represent as an unsigned") );
	unsigned val = 0;
	for (size_t i=0; i<value.size(); i++)
		{
		if (value[i] == true)
			{
			unsigned mask = 1 << i ;
			val |= mask;
			}
		}
	return val;
}


/** Initialize the discrete character state variable. We just set the value. */
void CharacterStateDiscrete::initialize(const std::vector<RbObject*> &attributes) {
    
    // set the state
    setState(static_cast<RbString*>( attributes[0] )->getValue()[0]);
//    numStates = static_cast<Natural*>( (RbLanguageObject*)(*attr)[0] )->getValue();
}

/** Is the character missing or ambiguous */
bool CharacterStateDiscrete::isMissingOrAmbiguous(void) const {

    if ( getNumberObservedStates() > 1 )
        return true;
    return false;
}


/** Set the number of states */
void CharacterStateDiscrete::setNumberOfStates(size_t n) {

    numStates = n;
    value.resize(n);
}


/** Set the values */
void CharacterStateDiscrete::setValue(const std::vector<bool> &v) {

    value = v;
}


/** Set value with single index */
void CharacterStateDiscrete::setValue(int x) {

    // reset all values
    for (size_t i=0; i<value.size(); i++) 
        {
        value[i] = false;
        }
    
    // set the value with given index to true
    value[x] = true;
}

