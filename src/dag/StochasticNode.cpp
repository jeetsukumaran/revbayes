/**
 * @file
 * This file contains the implementation of StochasticNode, which is derived
 * from DAGNode. StochasticNode is used for DAG nodes holding stochastic
 * variables with an associated probability density function.
 *
 * @brief Implementation of StochasticNode
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes development core team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-08-16, version 1.0
 * @extends DAGNode
 *
 * $Id$
 */

#include "Distribution.h"
#include "RbException.h"
#include "RbMove.h"
#include "RbMoveSchedule.h"
#include "StochasticNode.h"


/**
 * @brief Constructor from distribution
 *
 * The constructor sets the distribution and then draws a
 * random value, which is used as the starting value.
 *
 * @param dist      The distribution
 *
 */
StochasticNode::StochasticNode(Distribution* d, RandomNumberGenerator* r)
    : DAGNode(), clamped(false) {
    distribution = d;
    rng = r;
    value = distribution->rv(rng);
    storedValue = value->clone();
}   


/**
 * @brief Copy constructor
 *
 * Copy constructor for stochastic nodes.
 *
 * @param s     The stochastic node to be copied
 *
 */
StochasticNode::StochasticNode(const StochasticNode& s)
    : DAGNode(s), clamped(s.clamped), distribution(s.distribution) {
}


/**
 * @brief Destructor
 *
 * The stochastic node manages the distribution; the value is
 * managed and destroyed by the DAGNode base class.
 *
 */
StochasticNode::~StochasticNode() {

	if ( distribution != NULL )
		delete distribution;
}

RbObject& StochasticNode::operator=(const RbObject& obj) {

    try {
        // Use built-in fast down-casting first
        const StochasticNode& x = dynamic_cast<const StochasticNode&> (obj);

        StochasticNode& y = (*this);
        y = x;
        return y;
    } catch (std::bad_cast & bce) {
        try {
            // Try converting the value to an argumentRule
            const StochasticNode& x = dynamic_cast<const StochasticNode&> (*(obj.convertTo("const_node")));

            StochasticNode& y = (*this);
            y = x;
            return y;
        } catch (std::bad_cast & bce) {
            RbException e("Not supported assignment of " + obj.getClass()[0] + " to const_node");
            throw e;
        }
    }

    // dummy return
    return (*this);
}

DAGNode& StochasticNode::operator=(const DAGNode& obj) {

    try {
        // Use built-in fast down-casting first
        const StochasticNode& x = dynamic_cast<const StochasticNode&> (obj);

        StochasticNode& y = (*this);
        y = x;
        return y;
    } catch (std::bad_cast & bce) {
        try {
            // Try converting the value to an argumentRule
            const StochasticNode& x = dynamic_cast<const StochasticNode&> (*(obj.convertTo("const_node")));

            StochasticNode& y = (*this);
            y = x;
            return y;
        } catch (std::bad_cast & bce) {
            RbException e("Not supported assignment of " + obj.getClass()[0] + " to const_node");
            throw e;
        }
    }

    // dummy return
    return (*this);
}

StochasticNode& StochasticNode::operator=(const StochasticNode& obj) {

    changed = obj.changed;
    children = obj.children;
    (*lastMove) = (*obj.lastMove);
    monitors = obj.monitors;
    (*moves) = (*obj.moves);
    parents = obj.parents;
    (*storedValue) = (*obj.storedValue);
    touchedProbability = obj.touchedProbability;
    touchedLikelihood = obj.touchedLikelihood;
    (*value) = (*obj.value);

    return (*this);
}

/**
 * @brief Thouch affected nodes
 *
 * This function touches all affected DAG nodes, i.e. marks them as changed.
 *
 */
void StochasticNode::touchAffectedChildren() {

}

/**
 * @brief Thouch affected nodes
 *
 * This function touches all affected DAG nodes, i.e. marks them as changed.
 *
 */
void StochasticNode::touchAffectedParents() {

}

/**
 * @brief Tell affected DAG nodes to keep current value
 *
 * This function calls all affected DAG nodes so that they
 * have a chance to keep the current value and discard the
 * previous value.
 *
 */
void StochasticNode::keepAffectedChildren() {

}

/**
 * @brief Tell affected DAG nodes to keep current value
 *
 * This function calls all affected DAG nodes so that they
 * have a chance to keep the current value and discard the
 * previous value.
 *
 */
void StochasticNode::keepAffectedParents() {

}

/**
 * @brief Restore affected nodes
 *
 * This function calls all nodes that are affected by this DAG node and restores
 * them.
 *
 */
void StochasticNode::restoreAffectedChildren() {

}

/**
 * @brief Restore affected nodes
 *
 * This function calls all nodes that are affected by this DAG node and restores
 * them.
 *
 */
void StochasticNode::restoreAffectedParents() {

}

std::string StochasticNode::toString() const {
    return "Stochastic Node: " + value->toString();
}

/**
 * @brief Clamp the node to a value
 *
 * This function sets the value of the stochastic node to some
 * observed value (data). Then the clamped flag is set to make
 * sure that the node is marked as clamped and no change is made
 * to the observed value.
 *
 * @param observedVal   The observed value
 */
void StochasticNode::clamp(RbObject* observedVal) {

	if (storedValue != NULL)
        delete storedValue;
	if (value != NULL)
        delete value;
    storedValue = NULL;
    value = observedVal;
    clamped = true;
}

double StochasticNode::getLnProbabilityRatio(void) {

	double lnNumerator = distribution->lnPdf(this->value);
	double lnDenominator = distribution->lnPdf(this->storedValue);
	return lnNumerator - lnDenominator;
}

double StochasticNode::getLnProbability(void) {

	if (touchedProbability) {
		double lnProb = distribution->lnPdf(value);
		currentProbability = lnProb;
		touchedProbability = false;
	}
	return currentProbability;
}


void StochasticNode::initializeValue(RbObject* v) {

	value = v;
}

