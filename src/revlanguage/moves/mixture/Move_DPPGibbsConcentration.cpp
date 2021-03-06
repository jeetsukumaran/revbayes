#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DPPScaleCatValsMove.h"
#include "Move_DPPGibbsConcentration.h"
#include "DPPGibbsConcentrationMove.h"
#include "Natural.h"
#include "RbException.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"


using namespace RevLanguage;

Move_DPPGibbsConcentration::Move_DPPGibbsConcentration() : Move() {
    
}


/** Clone object */
Move_DPPGibbsConcentration* Move_DPPGibbsConcentration::clone(void) const {
    
	return new Move_DPPGibbsConcentration(*this);
}


void Move_DPPGibbsConcentration::constructInternalObject( void ) {
    // we free the memory first
    delete value;
    
    // now allocate a new vector-scale move
    double ne = static_cast<const RealPos &>( numElements->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<double>* tmp = static_cast<const RealPos &>( cp->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode< double > *sn = static_cast<RevBayesCore::StochasticNode<double> *>( tmp );
    RevBayesCore::TypedDagNode<int>* tmpNC = static_cast<const Integer &>( numCats->getRevObject() ).getDagNode();
    RevBayesCore::DeterministicNode< int > *nc = static_cast<RevBayesCore::DeterministicNode<int> *>( tmpNC );
    RevBayesCore::TypedDagNode<double>* gS = static_cast<const RealPos &>( gammaShape->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* gR = static_cast<const RealPos &>( gammaRate->getRevObject() ).getDagNode();

    
    value = new RevBayesCore::DPPGibbsConcentrationMove(sn, nc, gS, gR, ne, w);
}


/** Get class name of object */
const std::string& Move_DPPGibbsConcentration::getClassType(void) { 
    
    static std::string revClassType = "Move_DPPGibbsConcentration";
    
	return revClassType; 
}

/** Get class type spec describing type of object */
const TypeSpec& Move_DPPGibbsConcentration::getClassTypeSpec(void) { 
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return revClassTypeSpec; 
}



/** Return member rules (no members) */
const MemberRules& Move_DPPGibbsConcentration::getParameterRules(void) const {
    
    static MemberRules dppMove;
    static bool rulesSet = false;
    
    if ( !rulesSet )
    {
        
        dppMove.push_back( new ArgumentRule( "concentration", RealPos::getClassTypeSpec(), ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        dppMove.push_back( new ArgumentRule( "numDPPCats"   , Integer::getClassTypeSpec(), ArgumentRule::BY_CONSTANT_REFERENCE ) );
        dppMove.push_back( new ArgumentRule( "gammaShape"   , RealPos::getClassTypeSpec(), ArgumentRule::BY_CONSTANT_REFERENCE ) );
        dppMove.push_back( new ArgumentRule( "gammaRate"    , RealPos::getClassTypeSpec(), ArgumentRule::BY_CONSTANT_REFERENCE ) );
        dppMove.push_back( new ArgumentRule( "numElements"  , RealPos::getClassTypeSpec(), ArgumentRule::BY_CONSTANT_REFERENCE ) );
       
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        dppMove.insert( dppMove.end(), inheritedRules.begin(), inheritedRules.end() ); 
        
        rulesSet = true;
    }
    
    return dppMove;
}

/** Get type spec */
const TypeSpec& Move_DPPGibbsConcentration::getTypeSpec( void ) const {
    
    static TypeSpec typeSpec = getClassTypeSpec();
    
    return typeSpec;
}


/** Get type spec */
void Move_DPPGibbsConcentration::printValue(std::ostream &o) const {
    
    o << "Move_DPPGibbsConcentration(";
    if (cp != NULL) {
        o << cp->getName();
    }
    else if (numCats != NULL) {
        o << numCats->getName();
    }
    else if (gammaShape != NULL) {
        o << gammaShape->getName();
    }
    else if (gammaRate != NULL) {
        o << gammaRate->getName();
    }
    else if (numElements != NULL) {
        o << numElements->getName();
    }
    else {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Move_DPPGibbsConcentration::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "concentration" ) {
        cp = var;
    }
    else if ( name == "numDPPCats" ) {
        numCats = var;
    }
    else if ( name == "gammaShape" ) {
        gammaShape = var;
    }
    else if ( name == "gammaRate" ) {
        gammaRate = var;
    }
    else if ( name == "numElements" ) {
        numElements = var;
    }
    else {
        Move::setConstParameter(name, var);
    }
}
