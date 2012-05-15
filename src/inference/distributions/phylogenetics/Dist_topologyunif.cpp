///**
// * @file
// * This file contains the implementation of Dist_topologyunif, which holds
// * the parameters and functions related to a uniform distribution on
// * topologies (unrooted or rooted, binary or polytomous).
// *
// * @brief Implementation of Dist_topologyunif
// *
// * (c) Copyright 2009- under GPL version 3
// * @date Last modified: $Date$
// * @author The RevBayes Development Core Team
// * @license GPL version 3
// * @version 1.0
// * @since 2009-08-27, version 1.0
// *
// * $Id$
// */
//
//#include "RbBoolean.h"
//#include "DAGNode.h"
//#include "Dist_topologyunif.h"
//#include "Integer.h"
//#include "Natural.h"
//#include "RandomNumberFactory.h"
//#include "RandomNumberGenerator.h"
//#include "RbConstants.h"
//#include "RbNullObject.h"
//#include "RbException.h"
//#include "RbString.h"
//#include "RbUtil.h"
//#include "RbVector.h"
//#include "Simplex.h"
//#include "Topology.h"
//#include "TopologyNode.h"
//#include "ValueRule.h"
//
//#include <cmath>
//#include <sstream>
//#include <vector>
//
//
///** Default constructor for parser use */
//Dist_topologyunif::Dist_topologyunif( void ) : DistributionDiscrete( getMemberRules() ),
//                                                numTaxa( NULL ),
//                                                isRooted( NULL ),
//                                                isBinary( NULL ),
//                                                tipNames( NULL ) {
//
//    // Precalculate probability of topology
//    // \TODO: this will crash because we haven't set the parameter of the distribution
////    calculateTopologyProb();
//}
//
//
//void Dist_topologyunif::buildRandomBinaryTree(std::vector<TopologyNode*> &tips, size_t numTaxa) {
//    
//    if (tips.size() < numTaxa) {
//        // Get the rng
//        RandomNumberGenerator* rng = GLOBAL_RNG;
//        
//        // randomly draw one node from the list of tips
//        size_t index = static_cast<size_t>( floor(rng->uniform01()*tips.size()) );
//        
//        // get the node from the list
//        TopologyNode* parent = tips.at(index);
//        
//        // remove the randomly drawn node from the list
//        tips.erase(tips.begin()+index);
//        
//        // add a left child
//        TopologyNode* leftChild = new TopologyNode(0);
//        parent->addChild(leftChild);
//        leftChild->setParent(parent);
//        tips.push_back(leftChild);
//        
//        // add a right child
//        TopologyNode* rightChild = new TopologyNode(0);
//        parent->addChild(rightChild);
//        rightChild->setParent(parent);
//        tips.push_back(rightChild);
//        
//        // recursive call to this function
//        buildRandomBinaryTree(tips, numTaxa);
//    }
//}
//
//
//
///** Precalculate the probability and log probability of a topology drawn from the distribution. */
//void Dist_topologyunif::calculateTopologyProb( void ) {
//    
//    // Get the parameters
//    int  nTaxa  = static_cast<const Natural&  >( numTaxa->getValue()  ).getValue();
//    bool rooted = static_cast<const RbBoolean&>( isRooted->getValue() ).getValue();
//    bool binary = static_cast<const RbBoolean&>( isBinary->getValue() ).getValue();
//    
//    // Calculate probabilities of a topology drawn from the distribution
//    if (binary) {
//        // we have a binary topology
//        if (rooted) {
//            // set the topology probability to default values
//            lnTopologyProb = 0.0;
//            topologyProb = 1.0;
//            for (int i=2; i<=nTaxa; i++) {
//                double tmp = 1.0/(2.0*i - 3.0);
//                lnTopologyProb += log(tmp);
//                topologyProb *= tmp;
//            }
//        } else {
//            
//            // set the topology probability to default values
//            lnTopologyProb = 0.0;
//            topologyProb = 1.0;
//            for (int i=3; i<=nTaxa; i++) {
//                double tmp = 1.0/(2.0*i - 5.0);
//                lnTopologyProb += log(tmp);
//                topologyProb *= tmp;
//            }
//        }
//    }
//}
//
///** Precalculate the probability and log probability of a topology drawn from the distribution. */
//void Dist_topologyunif::calculateNumberOfStates( void ) {
//    
//    // Get the parameters
//    int  nTaxa  = static_cast<const Natural&  >( numTaxa->getValue()  ).getValue();
//    bool rooted = static_cast<const RbBoolean&>( isRooted->getValue() ).getValue();
//    bool binary = static_cast<const RbBoolean&>( isBinary->getValue() ).getValue();
//    
//    // Calculate probabilities of a topology drawn from the distribution
//    if (binary) {
//        // we have a binary topology
//        if (rooted) {
//            // set the number of states to default values
//            numberOfStates = 1;
//            for (int i=2; i<=nTaxa; i++) {
//                size_t tmp = (2*i - 3);
//                if (RbConstants::Size_t::max/numberOfStates < tmp) {
//                    numberOfStates = RbConstants::Size_t::inf;
//                    break;
//                }
//                numberOfStates *= tmp;
//            }
//        } else {
//            // set the number of states to default values
//            numberOfStates = 1;
//            for (int i=3; i<=nTaxa; i++) {
//                size_t tmp = (2*i - 5);
//                if (RbConstants::Size_t::max/numberOfStates < tmp) {
//                    numberOfStates = RbConstants::Size_t::inf;
//                    break;
//                }
//                numberOfStates *= tmp;
//            }
//        }
//    }
//}
//
//
///** Clone this object */
//Dist_topologyunif* Dist_topologyunif::clone( void ) const {
//
//    return new Dist_topologyunif( *this );
//}
//
//
///** Get class name of object */
//const std::string& Dist_topologyunif::getClassName(void) { 
//    
//    static std::string rbClassName = "Uniform topology distribution";
//    
//	return rbClassName; 
//}
//
///** Get class type spec describing type of object */
//const TypeSpec& Dist_topologyunif::getClassTypeSpec(void) { 
//    
//    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( DistributionDiscrete::getClassTypeSpec() ) );
//    
//	return rbClass; 
//}
//
///** Get type spec */
//const TypeSpec& Dist_topologyunif::getTypeSpec( void ) const {
//    
//    static TypeSpec typeSpec = getClassTypeSpec();
//    
//    return typeSpec;
//}
//
//
///** Get member variable rules */
//const MemberRules& Dist_topologyunif::getMemberRules( void ) const {
//
//    static MemberRules memberRules = MemberRules();
//    static bool        rulesSet = false;
//
//    if ( !rulesSet )
//		{
//        memberRules.push_back( new ValueRule( "numberTaxa" , Natural::getClassTypeSpec()      ) );
//        memberRules.push_back( new ValueRule( "tipNames"   , RbVector<RbString>::getClassTypeSpec() ) );
//        memberRules.push_back( new ValueRule( "isRooted"   , RbBoolean::getClassTypeSpec()    ) );
//        memberRules.push_back( new ValueRule( "isBinary"   , new RbBoolean(true) ) );
//
//        rulesSet = true;
//		}
//
//    return memberRules;
//}
//
//
///** Get the number of states in the distribution */
//size_t Dist_topologyunif::getNumberOfStates( void ) const {
//
//    return numberOfStates;
//}
//
//
///** Get the probability mass vector */
//const Simplex& Dist_topologyunif::getProbabilityMassVector( void ) {
//
//    int nTaxa  = static_cast<const Natural&>( numTaxa->getValue() ).getValue();
//
//    if ( nTaxa <= 6 ) {
//        probMassVector = Simplex( getNumberOfStates() );
//        return probMassVector;
//    }
//    else
//        throw RbException( "Too many topologies to generate probability mass vector" );
//}
//
//
///** Get the state vector for this distribution */
//const std::vector<RbLanguageObject*>& Dist_topologyunif::getStateVector( void ) const {
//    throw RbException("getStateVector not implemented for uniform topology distribution.");
//}
//
//
///** Get random variable type */
//const TypeSpec& Dist_topologyunif::getVariableType( void ) const {
//
//    static TypeSpec varTypeSpec = Topology::getClassTypeSpec();
//    
//    return varTypeSpec;
//}
//
//
///**
// * This function calculates the natural log of the
// * probability of a random topology variable.
// *
// * @brief Natural log of probability of topology
// *
// * @param value Observed topology
// * @return      Natural log of the probability
// */
//double Dist_topologyunif::lnPdf( const RbLanguageObject& value ) const {
//
//    // Get the parameters
//    const Topology& top = static_cast<const Topology&>( value );
//    int  nTaxa  = static_cast<const Natural&  >( numTaxa->getValue()  ).getValue();
//    bool rooted = static_cast<const RbBoolean&>( isRooted->getValue() ).getValue();
//    bool binary = static_cast<const RbBoolean&>( isBinary->getValue() ).getValue();
//
//    // If numTaxa, isRooted and isBinary fits top, then return precalculated probability stored in the class
//    // Otherwise return negative infinity
//    if ( top.getNumberOfNodes() == nTaxa && top.getIsRooted() == rooted && top.getIsBinary() == binary )
//        return lnTopologyProb;
//    else
//        return RbConstants::Double::neginf;
//}
//
//
///**
// * This function calculates the probability
// * of a categorical random variable.
// *
// * @brief Categorical probability density
// *
// * @param value Observed value
// * @return      Probability density
// */
//double Dist_topologyunif::pdf( const RbLanguageObject& value ) const {
//
//    // Get the parameters
//    const Topology& top = static_cast<const Topology&>( value );
//    int  nTaxa  = static_cast<const Natural&  >( numTaxa->getValue()  ).getValue();
//    bool rooted = static_cast<const RbBoolean&>( isRooted->getValue() ).getValue();
//    bool binary = static_cast<const RbBoolean&>( isBinary->getValue() ).getValue();
//
//    // If numTaxa, isRooted and isBinary fits top, then return precalculated probability stored in the class
//    // Otherwise return 0.0
//    if ( top.getNumberOfNodes() == nTaxa && top.getIsRooted() == rooted && top.getIsBinary() == binary )
//        return topologyProb;
//    else
//        return 0.0;
//}
//
//
///**
// * This function generates a topology from the distribution
// *
// * @brief Random draw from distribution
// *
// * @return      Randomly drawn topology
// */
//const RbLanguageObject& Dist_topologyunif::rv( void ) {
//
//    // Get the rng
//    RandomNumberGenerator* rng = GLOBAL_RNG;
//    
//    // Get the parameters
//    int  nTaxa                      = static_cast<const Natural&            >( numTaxa->getValue()  ).getValue();
//    bool rooted                     = static_cast<const RbBoolean&          >( isRooted->getValue() ).getValue();
//    bool binary                     = static_cast<const RbBoolean&          >( isBinary->getValue() ).getValue();
//    const RbVector<RbString>& names = static_cast<const RbVector<RbString>& >( tipNames->getValue() );
//
//    // Draw a random topology
//    if (binary) {
//        // internally we treat unrooted topologies the same as rooted
//        randomVariable.setIsRooted(rooted);
//        
//        TopologyNode* root = new TopologyNode((int)pow(2.0,nTaxa)-1);
//        std::vector<TopologyNode* > nodes;
//        nodes.push_back(root);
//        // recursively build the tree
//        buildRandomBinaryTree(nodes, nTaxa);
//        
//        // set tip names
//        for (int i=0; i<nTaxa; i++) {
//            size_t index = size_t( floor(rng->uniform01() * nodes.size()) );
//            
//            // get the node from the list
//            TopologyNode* node = nodes.at(index);
//            
//            // remove the randomly drawn node from the list
//            nodes.erase(nodes.begin()+index);
//            
//            // set name
//            std::string name = names[i];
//            node->setName(name);
//        }
//        
//        // initialize the topology by setting the root
//        randomVariable.setRoot(root);
//        
//        return randomVariable;
//    }
//    // \TODO: Draw a random multifurcating topology
//    return RbNullObject::getInstance();
//}
//
//
///** We intercept a call to set a member variable to make sure that the topology prob gets recalculated */
//void Dist_topologyunif::setMemberVariable( const std::string& name, const Variable* var ) {
//
//    if ( name == "numberTaxa" ) {
//        numTaxa = var;
//    } 
//    else if ( name == "isRooted" ) {
//        isRooted = var;
//    }
//    else if ( name == "isBinary" ) {
//        isBinary = var;
//    }
//    else if ( name == "tipNames" ) {
//        tipNames = var;
//    }
//    else {
//        DistributionDiscrete::setMemberVariable( name, var );
//    }
//    
//}
//