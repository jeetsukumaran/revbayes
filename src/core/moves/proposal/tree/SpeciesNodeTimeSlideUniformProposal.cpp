#include "SpeciesNodeTimeSlideUniformProposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "TreeUtilities.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
SpeciesNodeTimeSlideUniformProposal::SpeciesNodeTimeSlideUniformProposal( StochasticNode<TimeTree> *sp, std::vector< StochasticNode<TimeTree> *> gt ) : Proposal(),
    speciesTree( sp ),
    geneTrees( gt )
{
    // tell the base class to add the node
    addNode( speciesTree );
    
    for (size_t i=0; i < geneTrees.size(); ++i)
    {
        addNode( geneTrees[i] );
    }
    
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void SpeciesNodeTimeSlideUniformProposal::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
SpeciesNodeTimeSlideUniformProposal* SpeciesNodeTimeSlideUniformProposal::clone( void ) const
{
    
    return new SpeciesNodeTimeSlideUniformProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& SpeciesNodeTimeSlideUniformProposal::getProposalName( void ) const
{
    static std::string name = "SpeciesNodeTimeSlideUniform";
    
    return name;
}


/**
 * Perform the proposal.
 *
 * \return The hastings ratio.
 */
double SpeciesNodeTimeSlideUniformProposal::doProposal( void )
{
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    TimeTree& tau = speciesTree->getValue();
    
    // pick a random node which is not the root and neither the direct descendant of the root
    TopologyNode* node;
    do {
        double u = rng->uniform01();
        size_t index = size_t( std::floor(tau.getNumberOfNodes() * u) );
        node = &tau.getNode(index);
    } while ( node->isRoot() || node->isTip() );
    
    TopologyNode& parent = node->getParent();
    
    // we need to work with the times
    double parent_age  = parent.getAge();
    double my_age      = node->getAge();
    double child_Age   = node->getChild( 0 ).getAge();
    if ( child_Age < node->getChild( 1 ).getAge()) {
        child_Age = node->getChild( 1 ).getAge();
    }
    
    // now we store all necessary values
    storedNode = node;
    storedAge = my_age;
    
    // draw new ages and compute the hastings ratio at the same time
    double my_new_age = (parent_age-child_Age) * rng->uniform01() + child_Age;
    
    int upslideNodes = 0;
    int downslideNodes = 0;
    
    for ( size_t i=0; i<geneTrees.size(); ++i )
    {
        // get the i-th gene tree
        TimeTree& geneTree = geneTrees[i]->getValue();
        
        std::vector<TopologyNode*> nodes = getNodesInPopulation(geneTree, *node );
        
        for (size_t j=0; j<nodes.size(); ++j)
        {
            
            double a = nodes[j]->getAge();
            double new_a = a;
            if ( a > my_age )
            {
                ++upslideNodes;
                new_a = parent_age - (parent_age - my_new_age)/(parent_age - my_age) * (parent_age - a);
            }
            else
            {
                ++downslideNodes;
                new_a = child_Age + (my_new_age - child_Age)/(my_age - child_Age) * (a - child_Age);
            }
            
            // set the new age of this gene tree node
            geneTree.setAge( nodes[j]->getIndex(), new_a );
        }
        
    }
    
    // Sebastian: We need to work on a mechanism to make these proposal safe for non-ultrametric trees!
    //    if (min_age != 0.0)
    //    {
    //        for (size_t i = 0; i < tau.getNumberOfTips(); i++)
    //        {
    //            if (tau.getNode(i).getAge() < 0.0)
    //            {
    //                return RbConstants::Double::neginf;
    //            }
    //        }
    //    }
    
    
    // set the age of the species tree node
    tau.setAge( node->getIndex(), my_new_age );
    
    // compute the Hastings ratio
    double lnHastingsratio = upslideNodes * log( (parent_age - my_new_age)/(parent_age - my_age) ) + downslideNodes * log( (my_new_age - child_Age)/(my_age - child_Age) );
    
    return lnHastingsratio;
    
}


std::vector<TopologyNode*> SpeciesNodeTimeSlideUniformProposal::getNodesInPopulation( TimeTree &tau, TopologyNode &n )
{
    
    // I need all the oldest nodes/subtrees that have the same tips.
    // Those nodes need to be scaled too.
    
    // get the beginning and ending age of the population
    double max_age = -1.0;
    if ( n.isRoot() == false )
    {
        max_age = n.getParent().getAge();
    }
    
    // get all the taxa from the species tree that are descendants of node i
    double min_age_left = n.getChild(0).getAge();
    std::vector<TopologyNode*> speciesTaxa_left;
    TreeUtilities::getTaxaInSubtree( &n.getChild(0), speciesTaxa_left );
    
    // get all the individuals
    std::set<TopologyNode*> individualTaxa_left;
    for (size_t i = 0; i < speciesTaxa_left.size(); ++i)
    {
        const std::string &name = speciesTaxa_left[i]->getName();
        std::vector<TopologyNode*> ind = tau.getTipNodesWithSpeciesName( name );
        for (size_t j = 0; j < ind.size(); ++j)
        {
            individualTaxa_left.insert( ind[j] );
        }
    }
    
    // create the set of the nodes within this population
    std::set<TopologyNode*> nodesInPopulationSet;
    
    // now go through all nodes in the gene
    while ( individualTaxa_left.empty() == false )
    {
        std::set<TopologyNode*>::iterator it = individualTaxa_left.begin();
        individualTaxa_left.erase( it );
        
        TopologyNode *geneNode = *it;
        
        // add this node to our list of node we need to scale, if:
        // a) this is the root node
        // b) this is not the root and the age of the parent node is larger than the parent's age of the species node
        if ( geneNode->getAge() > min_age_left && geneNode->getAge() < max_age && geneNode->isTip() == false )
        {
            // add this node if it is within the age of our population
            nodesInPopulationSet.insert( geneNode );
        }
        
        if ( geneNode->isRoot() == false && ( max_age == -1.0 || max_age > geneNode->getParent().getAge() ) )
        {
            // push the parent to our current list
            individualTaxa_left.insert( &geneNode->getParent() );
        }
        
    }
    
    // get all the taxa from the species tree that are descendants of node i
    double min_age_right = n.getChild(1).getAge();
    std::vector<TopologyNode*> speciesTaxa_right;
    TreeUtilities::getTaxaInSubtree( &n.getChild(1), speciesTaxa_right );
    
    // get all the individuals
    std::set<TopologyNode*> individualTaxa_right;
    for (size_t i = 0; i < speciesTaxa_right.size(); ++i)
    {
        const std::string &name = speciesTaxa_right[i]->getName();
        std::vector<TopologyNode*> ind = tau.getTipNodesWithSpeciesName( name );
        for (size_t j = 0; j < ind.size(); ++j)
        {
            individualTaxa_right.insert( ind[j] );
        }
    }
    
    // now go through all nodes in the gene
    while ( individualTaxa_right.empty() == false )
    {
        std::set<TopologyNode*>::iterator it = individualTaxa_right.begin();
        individualTaxa_right.erase( it );
        
        TopologyNode *geneNode = *it;
        
        // add this node to our list of node we need to scale, if:
        // a) this is the root node
        // b) this is not the root and the age of the parent node is larger than the parent's age of the species node
        if ( geneNode->getAge() > min_age_right && geneNode->getAge() < max_age && geneNode->isTip() == false )
        {
            // add this node if it is within the age of our population
            nodesInPopulationSet.insert( geneNode );
        }
        
        if ( geneNode->isRoot() == false && ( max_age == -1.0 || max_age > geneNode->getParent().getAge() ) )
        {
            // push the parent to our current list
            individualTaxa_right.insert( &geneNode->getParent() );
        }
        
    }
    
    
    
    
    // convert the set into a vector
    std::vector<TopologyNode*> nodesInPopulation;
    for (std::set<TopologyNode*>::iterator it = nodesInPopulationSet.begin(); it != nodesInPopulationSet.end(); ++it)
    {
        nodesInPopulation.push_back( *it );
    }
    
    return nodesInPopulation;
}


/**
 *
 */
void SpeciesNodeTimeSlideUniformProposal::prepareProposal( void )
{
    
}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void SpeciesNodeTimeSlideUniformProposal::printParameterSummary(std::ostream &o) const
{
    
    // no parameters
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void SpeciesNodeTimeSlideUniformProposal::undoProposal( void )
{
    
    // undo the proposal
    
    
    TopologyNode& parent = storedNode->getParent();
    
    // we need to work with the times
    double parent_age  = parent.getAge();
    double my_new_age      = storedNode->getAge();
    double child_Age   = storedNode->getChild( 0 ).getAge();
    if ( child_Age < storedNode->getChild( 1 ).getAge())
    {
        child_Age = storedNode->getChild( 1 ).getAge();
    }
    
    for ( size_t i=0; i<geneTrees.size(); ++i )
    {
        // get the i-th gene tree
        TimeTree& geneTree = geneTrees[i]->getValue();
        
        std::vector<TopologyNode*> nodes = getNodesInPopulation(geneTree, *storedNode );
        
        for (size_t j=0; j<nodes.size(); ++j)
        {
            
            double new_a = nodes[j]->getAge();
            double a = new_a;
            if ( new_a > my_new_age )
            {
                a = parent_age - (parent_age - storedAge)/(parent_age - my_new_age) * (parent_age - new_a);
            }
            else
            {
                a = child_Age + (storedAge - child_Age)/(my_new_age - child_Age) * (new_a - child_Age);
            }
            
            // set the new age of this gene tree node
            geneTree.setAge( nodes[j]->getIndex(), a );
        }
        
    }
    
    // set the age of the species tree node
    speciesTree->getValue().setAge( storedNode->getIndex(), storedAge );
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void SpeciesNodeTimeSlideUniformProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if ( oldN == speciesTree )
    {
        speciesTree = static_cast<StochasticNode<TimeTree>* >(newN) ;
    }
    else
    {
        for ( size_t i=0; i<geneTrees.size(); ++i )
        {
            if ( oldN == geneTrees[i] )
            {
                geneTrees[i] = static_cast<StochasticNode<TimeTree>* >(newN) ;
            }
        }
    }
    
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void SpeciesNodeTimeSlideUniformProposal::tune( double rate )
{
    
    // nothing to tune
    
}

