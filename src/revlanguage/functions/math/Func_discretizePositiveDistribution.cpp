#include "ArgumentRule.h"
#include "DiscretizeDistributionFunction.h"
#include "Ellipsis.h"
#include "Func_discretizePositiveDistribution.h"
#include "Integer.h"
#include "ModelVector.h"
#include "RbUtil.h"
#include "Real.h"
#include "RlBoolean.h"
#include "RlPositiveContinuousDistribution.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"


using namespace RevLanguage;

Func_discretizePositiveDistribution::Func_discretizePositiveDistribution() : TypedFunction< ModelVector< RealPos > >()
{
    
}

/* Clone object */
Func_discretizePositiveDistribution* Func_discretizePositiveDistribution::clone( void ) const
{
    
    return new Func_discretizePositiveDistribution( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >* Func_discretizePositiveDistribution::createFunction( void ) const
{
    
    const PositiveContinuousDistribution& rlDistribution = static_cast<const PositiveContinuousDistribution &>( this->args[0].getVariable()->getRevObject() );
    RevBayesCore::ContinuousDistribution* g0 = static_cast<RevBayesCore::ContinuousDistribution* >( rlDistribution.createDistribution() );
    RevBayesCore::TypedDagNode<int>* numCats = static_cast<const Integer &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::DiscretizeDistributionFunction *func = new RevBayesCore::DiscretizeDistributionFunction( g0, numCats );
    
    return func;
}


/** Get argument rules */
const ArgumentRules& Func_discretizePositiveDistribution::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rulesSet = false;
    
    if ( !rulesSet )
    {
        
        argumentRules.push_back( new ArgumentRule( "G0"      , PositiveContinuousDistribution::getClassTypeSpec(), ArgumentRule::BY_CONSTANT_REFERENCE ) );
        argumentRules.push_back( new ArgumentRule( "num_cats", Integer::getClassTypeSpec()               , ArgumentRule::BY_VALUE ) );
        rulesSet = true;
    }
    
    return argumentRules;
}


/** Get class name of object */
const std::string& Func_discretizePositiveDistribution::getClassName(void)
{
    
    static std::string rbClassName = "Func_discretizePositiveDistribution";
    
    return rbClassName;
}


/** Get class type spec describing type of object */
const RevLanguage::TypeSpec& Func_discretizePositiveDistribution::getClassTypeSpec(void) {
    
    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rbClass;
}


/** Get type spec */
const TypeSpec& Func_discretizePositiveDistribution::getTypeSpec( void ) const {
    
    static TypeSpec typeSpec = getClassTypeSpec();
    
    return typeSpec;
}
