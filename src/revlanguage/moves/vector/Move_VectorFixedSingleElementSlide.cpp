//
//  MoveSlide.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 8/6/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlBoolean.h"
#include "ContinuousStochasticNode.h"
#include "ModelVector.h"
#include "Move_VectorFixedSingleElementSlide.h"
#include "Natural.h"
#include "RbException.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "VectorFixedSingleElementSlidingMove.h"


using namespace RevLanguage;

Move_VectorFixedSingleElementSlide::Move_VectorFixedSingleElementSlide() : Move() {
    
}

/** Clone object */
Move_VectorFixedSingleElementSlide* Move_VectorFixedSingleElementSlide::clone(void) const {
    
	return new Move_VectorFixedSingleElementSlide(*this);
}


void Move_VectorFixedSingleElementSlide::constructInternalObject( void ) {

    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double l = static_cast<const RealPos &>( lambda->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* tmp = static_cast<const ModelVector<Real> &>( v->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> > *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> > *>( tmp );
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    size_t e = static_cast<const Natural &>( whichElement->getRevObject() ).getValue();
    value = new RevBayesCore::VectorFixedSingleElementSlidingMove(n, l, t, w, e-1);
}


/** Get Rev type of object */
const std::string& Move_VectorFixedSingleElementSlide::getClassType(void) { 
    
    static std::string revType = "Move_VectorFixedSingleElementSlide";
	return revType;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_VectorFixedSingleElementSlide::getClassTypeSpec(void) { 
    
    static TypeSpec revTypeSpec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
	return revTypeSpec; 
}



/** Return member rules (no members) */
const MemberRules& Move_VectorFixedSingleElementSlide::getParameterRules(void) const {
    
    static MemberRules moveMemberRules;
    static bool rulesSet = false;
    
    if ( !rulesSet )
        {
        moveMemberRules.push_back( new ArgumentRule( "x"      , ModelVector<Real>::getClassTypeSpec(), ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        moveMemberRules.push_back( new ArgumentRule( "lambda" , RealPos::getClassTypeSpec()          , ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new Real(1.0) ) );
        moveMemberRules.push_back( new ArgumentRule( "tune"   , RlBoolean::getClassTypeSpec()        , ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );
        moveMemberRules.push_back( new ArgumentRule( "element", Natural::getClassTypeSpec()          , ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new Natural( 1 ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        moveMemberRules.insert( moveMemberRules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rulesSet = true;
        }
    return moveMemberRules;
}

/** Get type spec */
const TypeSpec& Move_VectorFixedSingleElementSlide::getTypeSpec( void ) const {
    
    static TypeSpec typeSpec = getClassTypeSpec();
    
    return typeSpec;
}


/** Get type spec */
void Move_VectorFixedSingleElementSlide::printValue(std::ostream &o) const {
    
    o << "Move_VectorFixedSingleElementSlide(";
    if (v != NULL) {
        o << v->getName();
    }
    else {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Move_VectorFixedSingleElementSlide::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "x" ) {
        v = var;
    }
    else if ( name == "lambda" ) {
        lambda = var;
    }
    else if ( name == "weight" ) {
        weight = var;
    }
    else if ( name == "tune" ) {
        tune = var;
    }
    else if ( name == "element" ) {
        whichElement = var;
    }
    else {
        Move::setConstParameter(name, var);
    }
}
