#include "ConstantNode.h"
#include "Integer.h"
#include "Natural.h"
#include "RlBoolean.h"
#include "Probability.h"
#include "RlMatrixRealSymmetric.h"
#include "RbUtil.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "RlMemberFunction.h"

#include <iomanip>
#include <sstream>

using namespace RevLanguage;

/* Default constructor */
MatrixRealSymmetric::MatrixRealSymmetric(void) : MatrixReal() {
    
    ArgumentRules* precisionArgRules = new ArgumentRules();
    methods.addFunction("precision", new MemberProcedure( Natural::getClassTypeSpec(), precisionArgRules) );
}


/* Construct from double */
MatrixRealSymmetric::MatrixRealSymmetric( RevBayesCore::MatrixReal * mat ) : MatrixReal( mat ) {
    
    ArgumentRules* precisionArgRules = new ArgumentRules();
    methods.addFunction("precision", new MemberProcedure( Natural::getClassTypeSpec(), precisionArgRules) );
}


/* Construct from double */
MatrixRealSymmetric::MatrixRealSymmetric( RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal> * mat ) : MatrixReal( mat ) {
    
    ArgumentRules* precisionArgRules = new ArgumentRules();
    methods.addFunction("precision", new MemberProcedure( Natural::getClassTypeSpec(), precisionArgRules) );
}


/** Clone object */
MatrixRealSymmetric* MatrixRealSymmetric::clone(void) const {
    
	return new MatrixRealSymmetric(*this);
}


/** Convert to type. The caller manages the returned object. */
RevObject* MatrixRealSymmetric::convertTo( const TypeSpec& type ) const {
    
    return RevObject::convertTo( type );
}

RevPtr<RevVariable> MatrixRealSymmetric::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found) {
    
    if (name == "precision")
        {
        found = true;
        
        //int n = (int)this->dagNode->getValue().getNumberOfCharacters();

        RevBayesCore::MatrixReal& m = this->dagNode->getValue();
        std::cout << m[0][0] << std::endl;

        // get the member with give index
        /*const RevBayesCore::MatrixReal& m = static_cast<const RevBayesCore::MatrixReal&>( args[0].getVariable()->getRevObject() );
        
        if (this->dagNode->getValue().getNumberOfStates() < (size_t)(index.getValue()) ) {
            throw RbException("Index out of bounds in []");
            }
        
        const std::vector<double>& element = this->dagNode->getValue()[ size_t(index.getValue()) - 1];
        RevBayesCore::RbVector<double> elementVector;
        for (size_t i=0; i < this->dagNode->getValue().size(); ++i) {
            elementVector.push_back( element[i] );
            }
        
        return new RevVariable( new ModelVector<Real>( elementVector ) );*/
        }
    return ModelObject<RevBayesCore::MatrixReal>::executeMethod( name, args, found );
}

/** Get Rev type of object */
const std::string& MatrixRealSymmetric::getClassType(void) {
    
    static std::string revType = "MatrixRealSymmetric";
    
	return revType;
}

/** Get class type spec describing type of object */
const TypeSpec& MatrixRealSymmetric::getClassTypeSpec(void) {
    
    static TypeSpec revTypeSpec = TypeSpec( getClassType(), new TypeSpec( MatrixReal::getClassTypeSpec() ) );
    
	return revTypeSpec;
}

/** Get type spec */
const TypeSpec& MatrixRealSymmetric::getTypeSpec( void ) const {
    
    static TypeSpec typeSpec = getClassTypeSpec();
    
    return typeSpec;
}


/** Print value for user */
void MatrixRealSymmetric::printValue(std::ostream &o) const {
    
    long previousPrecision = o.precision();
    std::ios_base::fmtflags previousFlags = o.flags();
    
    std::fixed( o );
    o.precision( 3 );

    dagNode->printValue( o );
    
    o.setf( previousFlags );
    o.precision( previousPrecision );
}


