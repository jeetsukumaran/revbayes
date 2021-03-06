
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "GewekeStoppingRule.h"
#include "OptionRule.h"
#include "Probability.h"
#include "RlGewekeStoppingRule.h"
#include "RealPos.h"
#include "RbException.h"
#include "RlString.h"
#include "TypeSpec.h"


using namespace RevLanguage;

/**
 * Default constructor.
 * Create the default instance.
 */
GewekeStoppingRule::GewekeStoppingRule(void) : AbstractConvergenceStoppingRule()
{
    
}


/**
 * Clone object
 */
GewekeStoppingRule* GewekeStoppingRule::clone(void) const
{
    
    return new GewekeStoppingRule(*this);
}


void GewekeStoppingRule::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new stopping rule
    double p = static_cast<const Probability &>( prob->getRevObject() ).getValue();
    double f1 = static_cast<const Probability &>( frac1->getRevObject() ).getValue();
    double f2 = static_cast<const Probability &>( frac2->getRevObject() ).getValue();
    int fq = static_cast<const Natural &>( frequency->getRevObject() ).getValue();
    const std::string &fn = static_cast<const RlString &>( filename->getRevObject() ).getValue();
    
    RevBayesCore::BurninEstimatorContinuous *burninEst = constructBurninEstimator();
    
    value = new RevBayesCore::GewekeStoppingRule(p, f1, f2, fn, size_t(fq), burninEst);
    
}


/** Get Rev type of object */
const std::string& GewekeStoppingRule::getClassType(void)
{
    
    static std::string revType = "GewekeStoppingRule";
    
    return revType;
}

/** Get class type spec describing type of object */
const TypeSpec& GewekeStoppingRule::getClassTypeSpec(void)
{
    
    static TypeSpec revTypeSpec = TypeSpec( getClassType(), new TypeSpec( AbstractConvergenceStoppingRule::getClassTypeSpec() ) );
    
    return revTypeSpec;
}



/** Return member rules */
const MemberRules& GewekeStoppingRule::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rulesSet = false;
    
    if ( !rulesSet )
    {
        
        memberRules.push_back( new ArgumentRule( "prob" , Probability::getClassTypeSpec() , ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.05) ) );
        memberRules.push_back( new ArgumentRule( "frac1", Probability::getClassTypeSpec() , ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.1) ) );
        memberRules.push_back( new ArgumentRule( "frac2", Probability::getClassTypeSpec() , ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.5) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = AbstractConvergenceStoppingRule::getParameterRules();
        memberRules.insert( memberRules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        
        rulesSet = true;
    }
    
    return memberRules;
}

/** Get type spec */
const TypeSpec& GewekeStoppingRule::getTypeSpec( void ) const
{
    
    static TypeSpec typeSpec = getClassTypeSpec();
    
    return typeSpec;
}


/** Get type spec */
void GewekeStoppingRule::printValue(std::ostream &o) const
{
    
    o << "GewekeStoppingRule";
}


/** Set a member variable */
void GewekeStoppingRule::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "prob" )
    {
        prob = var;
    }
    else if ( name == "frac1" )
    {
        frac1 = var;
    }
    else if ( name == "frac2" )
    {
        frac2 = var;
    }
    else
    {
        AbstractConvergenceStoppingRule::setConstParameter(name, var);
    }
    
}
