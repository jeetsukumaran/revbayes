#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RbLanguageObject.h"
#include "RbException.h"
#include "RealPos.h"
#include "RlBoolean.h"
#include "RlRateAgeBetaShift.h"
#include "RlTimeTree.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "Vector.h"


using namespace RevLanguage;

RateAgeBetaShift::RateAgeBetaShift() : Move() {
    
}


/** Clone object */
RateAgeBetaShift* RateAgeBetaShift::clone(void) const {
    
	return new RateAgeBetaShift(*this);
}


void RateAgeBetaShift::constructInternalObject( void ) {
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    RevBayesCore::TypedDagNode<RevBayesCore::TimeTree> *tmp = static_cast<const TimeTree &>( tree->getValue() ).getValueNode();
    double d = static_cast<const RealPos &>( delta->getValue() ).getValue();
    bool at = static_cast<const RlBoolean &>( tune->getValue() ).getValue();
    double w = static_cast<const RealPos &>( weight->getValue() ).getValue();
    RevBayesCore::StochasticNode<RevBayesCore::TimeTree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::TimeTree> *>( tmp );
//    value = new RevBayesCore::RateAgeBetaShift(t, d, at, w);
}


/** Get class name of object */
const std::string& RateAgeBetaShift::getClassName(void) { 
    
    static std::string rbClassName = "RateAgeBetaShift";
    
	return rbClassName; 
}

/** Get class type spec describing type of object */
const TypeSpec& RateAgeBetaShift::getClassTypeSpec(void) { 
    
    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return rbClass; 
}



/** Return member rules (no members) */
const MemberRules& RateAgeBetaShift::getMemberRules(void) const {
    
    static MemberRules nniMemberRules;
    static bool rulesSet = false;
    
    if ( !rulesSet ) {
        nniMemberRules.push_back( new ArgumentRule( "tree", false, TimeTree::getClassTypeSpec() ) );
        nniMemberRules.push_back( new ArgumentRule( "rates", false, Vector<RealPos>::getClassTypeSpec() ) );
        nniMemberRules.push_back( new ArgumentRule( "delta", true, RealPos::getClassTypeSpec() , new Real(1.0) ) );
        nniMemberRules.push_back( new ArgumentRule( "tune"  , true, RlBoolean::getClassTypeSpec(), new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getMemberRules();
        nniMemberRules.insert( nniMemberRules.end(), inheritedRules.begin(), inheritedRules.end() ); 
        
        rulesSet = true;
    }
    
    return nniMemberRules;
}

/** Get type spec */
const TypeSpec& RateAgeBetaShift::getTypeSpec( void ) const {
    
    static TypeSpec typeSpec = getClassTypeSpec();
    
    return typeSpec;
}



/** Get type spec */
void RateAgeBetaShift::printValue(std::ostream &o) const {
    
    o << "RateAgeBetaShift(";
    if (tree != NULL) {
        o << tree->getName();
    }
    else {
        o << "?";
    }
    o << ")";
}


/** Set a NearestNeighborInterchange variable */
void RateAgeBetaShift::setConstMemberVariable(const std::string& name, const RbPtr<const Variable> &var) {
    
    if ( name == "tree" ) {
        tree = var;
    }
    else if ( name == "rates" ) {
        rates = var;
    }
    else if ( name == "delta" ) {
        delta = var;
    }
    else if ( name == "tune" ) {
        tune = var;
    }
    else {
        Move::setConstMemberVariable(name, var);
    }
}