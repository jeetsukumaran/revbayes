#include "TajimasDFunction.h"
#include "RbException.h"
#include "RbMathCombinatorialFunctions.h"

#include <cmath>

using namespace RevBayesCore;

TajimasDFunction::TajimasDFunction(const TypedDagNode<AbstractDiscreteCharacterData> *a) : TypedFunction<double>( new double(0.0) ),
    alignment( a )
{
    // add the lambda parameter as a parent
    addParameter( alignment );
    
    update();
}


TajimasDFunction::~TajimasDFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



TajimasDFunction* TajimasDFunction::clone( void ) const
{
    
    return new TajimasDFunction( *this );
}


void TajimasDFunction::update( void )
{
    int S = int( alignment->getValue().getNumberOfSegregatingSites() );
    size_t n = alignment->getValue().getNumberOfTaxa();
    
    double a1 = RbMath::harmonicNumber(n-1);
    double a2 = RbMath::squaredHarmonicNumber(n-1);
    
    double pi  = alignment->getValue().getPaiwiseSequenceDifference();
    double theta = S / a1;
    
    double b1 = (n+1.0)/ double(3.0*(n-1.0));
    double b2 = 2.0*(n*n+n+3.0) / double(9.0*n*(n-1.0));
    
    double c1 = b1 - 1.0/a1;
    double c2 = b2 - (n+2.0)/(a1*n) + a2/(a1*a1);

    double e1 = c1 / a1;
    double e2 = c2 / (a1*a1 + a2);
    
    double C = e1*S + e2*S*(S-1.0);
    
    *value = (pi - theta) / sqrt(C);
}



void TajimasDFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == alignment)
    {
        alignment = static_cast<const TypedDagNode< AbstractDiscreteCharacterData >* >( newP );
    }
    
}


