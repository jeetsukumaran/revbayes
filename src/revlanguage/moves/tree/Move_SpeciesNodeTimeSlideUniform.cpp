#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_SpeciesNodeTimeSlideUniform.h"
#include "SpeciesNodeTimeSlideUniformProposal.h"
#include "RbException.h"
#include "RealPos.h"
#include "RlTimeTree.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"


using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Move_SpeciesNodeTimeSlideUniform::Move_SpeciesNodeTimeSlideUniform() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the move.
 */
Move_SpeciesNodeTimeSlideUniform* Move_SpeciesNodeTimeSlideUniform::clone(void) const
{
    
    return new Move_SpeciesNodeTimeSlideUniform(*this);
}


/**
 * Create a new internal move object.
 *
 * This function simply dynamically allocates a new internal move object that is
 * associated with the variable (DAG-node). The internal move object is created by calling its
 * constructor and passing the move-parameters (the variable and other parameters) as arguments of the
 * constructor. The move constructor takes care of the proper hook-ups.
 *
 * \return A new internal distribution object.
 */
void Move_SpeciesNodeTimeSlideUniform::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::TimeTree>* tmp = static_cast<const TimeTree &>( speciesTree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::TimeTree> *st = static_cast<RevBayesCore::StochasticNode<RevBayesCore::TimeTree> *>( tmp );
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::TimeTree> >* tmp_gt = static_cast<const ModelVector<TimeTree> &>( geneTrees->getRevObject() ).getDagNode();
    std::set<const RevBayesCore::DagNode*> parents = tmp_gt->getParents();
    std::vector< RevBayesCore::StochasticNode<RevBayesCore::TimeTree> *> gt;
    for (std::set<const RevBayesCore::DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it)
    {
        const RevBayesCore::StochasticNode<RevBayesCore::TimeTree> *theNode = dynamic_cast< const RevBayesCore::StochasticNode<RevBayesCore::TimeTree>* >( *it );
        if ( theNode != NULL )
        {
            gt.push_back( const_cast< RevBayesCore::StochasticNode<RevBayesCore::TimeTree>* >( theNode ) );
        }
        else
        {
            throw RbException("Could not create a mvElementScale because the node isn't a vector of stochastic nodes.");
        }
    }
    
    
    RevBayesCore::Proposal *p = new RevBayesCore::SpeciesNodeTimeSlideUniformProposal(st,gt);
    value = new RevBayesCore::MetropolisHastingsMove(p,w);
    
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Move_SpeciesNodeTimeSlideUniform::getClassType(void)
{
    
    static std::string revType = "Move_SpeciesNodeTimeSlideUniform";
    
    return revType;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Move_SpeciesNodeTimeSlideUniform::getClassTypeSpec(void)
{
    
    static TypeSpec revTypeSpec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return revTypeSpec;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the scale move are:
 * (1) the variable which must be a time-tree.
 *
 * \return The member rules.
 */
const MemberRules& Move_SpeciesNodeTimeSlideUniform::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rulesSet = false;
    
    if ( !rulesSet )
    {
        memberRules.push_back( new ArgumentRule( "speciesTree", TimeTree::getClassTypeSpec()             , ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC    ) );
        memberRules.push_back( new ArgumentRule( "geneTrees"  , ModelVector<TimeTree>::getClassTypeSpec(), ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        memberRules.insert( memberRules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rulesSet = true;
    }
    
    return memberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Move_SpeciesNodeTimeSlideUniform::getTypeSpec( void ) const
{
    
    static TypeSpec typeSpec = getClassTypeSpec();
    
    return typeSpec;
}


/**
 * Print the value for the user.
 */
void Move_SpeciesNodeTimeSlideUniform::printValue(std::ostream &o) const
{
    
    o << "SpeciesNodeTimeSlideUniform(";
    if (speciesTree != NULL)
    {
        o << speciesTree->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
    
}


/**
 * Set a member variable.
 *
 * Sets a member variable with the given name and store the pointer to the variable.
 * The value of the variable might still change but this function needs to be called again if the pointer to
 * the variable changes. The current values will be used to create the distribution object.
 *
 * \param[in]    name     Name of the member variable.
 * \param[in]    var      Pointer to the variable.
 */
void Move_SpeciesNodeTimeSlideUniform::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "speciesTree" )
    {
        speciesTree = var;
    }
    else if ( name == "geneTrees" )
    {
        geneTrees = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
    
}




