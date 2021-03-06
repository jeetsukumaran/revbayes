/*
 * File:   PhyloOrnsteinUhlenbeckProcess.cpp
 * Author: Bastien Boussau
 *
 * Created on 04/10/14
 */

#include "PhyloOrnsteinUhlenbeckProcess.h"
#include "ConstantNode.h"
#include "DistributionNormal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "StochasticNode.h"

#include <cmath>
#include <iostream>



using namespace RevBayesCore;


// constructor(s)
PhyloOrnsteinUhlenbeckProcess::PhyloOrnsteinUhlenbeckProcess(const TypedDagNode< TimeTree > *t): TypedDistribution< RbVector<double> >(new RbVector<double>()),
    tau( t )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( tau );
    
    homogeneousSigma = NULL;
    homogeneousMean = NULL;
    homogeneousPhi = NULL;
    heterogeneousSigma = NULL;
    heterogeneousMean = NULL;
    heterogeneousPhi = NULL;
    root = NULL;
    
 //   simulate();
}



PhyloOrnsteinUhlenbeckProcess* PhyloOrnsteinUhlenbeckProcess::clone(void) const
{
    return new PhyloOrnsteinUhlenbeckProcess( *this );
}


double PhyloOrnsteinUhlenbeckProcess::computeLnProbability(void)
{
    
    double ret = recursiveLnProb(tau->getValue().getRoot());
    return ret;
}


double PhyloOrnsteinUhlenbeckProcess::recursiveLnProb( const TopologyNode& from ) {
    
    double lnProb = 0.0;
    size_t index = from.getIndex();
    double val = (*value)[index];
    
    if (from.isRoot())
    {

        // we only need to check that the root has the same value as set by the root parameter
        if ( val != root->getValue() )
        {
            return RbConstants::Double::neginf;
        }
        
    }
    else
    {
        
        // x ~ normal(x_up, sigma^2 * branchLength)
        double upval = (*value)[from.getParent().getIndex()];
        double t = from.getBranchLength();
        if ( homogeneousPhi != NULL ) // homogeneous
        {
            double e = exp(-homogeneousPhi->getValue() * t);
            double e2 = exp(-2 * homogeneousPhi->getValue() * t);
            double m = e * upval + (1 - e) * homogeneousMean->getValue();
            double standDev = homogeneousSigma->getValue() * sqrt((1 - e2) / 2 / homogeneousPhi->getValue());
        
            lnProb += RbStatistics::Normal::lnPdf(val, standDev, m);
        }
        else //heterogeneous
        {
            size_t nodeIndex = from.getIndex();
            double e = exp(-heterogeneousPhi->getValue()[nodeIndex] * t);
            double e2 = exp(-2 * heterogeneousPhi->getValue()[nodeIndex] * t);
            double m = e * upval + (1 - e) * heterogeneousMean->getValue()[nodeIndex];
            double standDev = heterogeneousSigma->getValue()[nodeIndex] * sqrt((1 - e2) / 2 / heterogeneousPhi->getValue()[nodeIndex]);
            
            lnProb += RbStatistics::Normal::lnPdf(val, standDev, m);

        }
    }
    
    // propagate forward
    size_t numChildren = from.getNumberOfChildren();
    
    for (size_t i = 0; i < numChildren; ++i)
    {
        lnProb += recursiveLnProb(from.getChild(i));
    }
    
    return lnProb;
    
}

void PhyloOrnsteinUhlenbeckProcess::redrawValue(void)
{
    simulate();
}


void PhyloOrnsteinUhlenbeckProcess::simulate()
{
    
    recursiveSimulate(tau->getValue().getRoot());
}


void PhyloOrnsteinUhlenbeckProcess::recursiveSimulate(const TopologyNode& from)  {
    
    size_t index = from.getIndex();
    
    if (from.isRoot())
    {

        // we can set directly the value of the root
        (*value)[index] = root->getValue();
    }
    else
    {
        
        // x ~ normal(x_up, sigma^2 * branchLength)
        
        double upval = (*value)[from.getParent().getIndex()];
        double t = from.getBranchLength();
        double m = 0.0;
        double standDev = 1.0;
        if ( homogeneousPhi != NULL ) // homogeneous
        {
            double e = exp(-homogeneousPhi->getValue() * t);
            double e2 = exp(-2 * homogeneousPhi->getValue() * t);
            m = e * upval + (1 - e) * homogeneousMean->getValue();
            standDev = homogeneousSigma->getValue() * sqrt((1 - e2) / 2 / homogeneousPhi->getValue());
        }
        else // heterogeneous
        {
            size_t nodeIndex = from.getIndex();
            double e = exp(-heterogeneousPhi->getValue()[nodeIndex] * t);
            double e2 = exp(-2 * heterogeneousPhi->getValue()[nodeIndex] * t);
            m = e * upval + (1 - e) * heterogeneousMean->getValue()[nodeIndex];
            standDev = heterogeneousSigma->getValue()[nodeIndex] * sqrt((1 - e2) / 2 / heterogeneousPhi->getValue()[nodeIndex]);
        }
        // simulate the new Val
        RandomNumberGenerator* rng = GLOBAL_RNG;
        (*value)[index] = RbStatistics::Normal::rv(m, standDev, *rng);
        
    }
    
    // propagate forward
    size_t numChildren = from.getNumberOfChildren();
    for (size_t i = 0; i < numChildren; ++i)
    {
        recursiveSimulate(from.getChild(i));
    }
    
}


/** Swap a parameter of the distribution */
void PhyloOrnsteinUhlenbeckProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if ( oldP == tau )
    {
        tau = static_cast< const TypedDagNode<TimeTree> * >( newP );
    }
    
    if ( oldP == homogeneousSigma )
    {
        homogeneousSigma = static_cast< const TypedDagNode<double> * >( newP );
    }
    
    if ( oldP == homogeneousMean ) {
        homogeneousMean = static_cast< const TypedDagNode< double > * >( newP );
    }
    
    if ( oldP == homogeneousPhi ) {
        homogeneousPhi = static_cast< const TypedDagNode< double > * >( newP );
    }
    
    if ( oldP == heterogeneousSigma ) {
        heterogeneousSigma = static_cast< const TypedDagNode< RbVector< double > > * >( newP );
    }
    
    if ( oldP == heterogeneousMean ) {
        heterogeneousMean = static_cast< const TypedDagNode< RbVector< double > > * >( newP );
    }
    
    if ( oldP == heterogeneousPhi ) {
        heterogeneousPhi = static_cast< const TypedDagNode< RbVector< double > > * >( newP );
    }
    
    if ( oldP == root ) {
        root = static_cast< const TypedDagNode<double> * >( newP );
    }
    
}


void PhyloOrnsteinUhlenbeckProcess::setMean( const TypedDagNode< double >* x)
{
    // remove the old parameter
    removeParameter( heterogeneousMean );
    removeParameter( homogeneousMean );
    
    // set our local pointer copy
    homogeneousMean = x;
    
    // add the new parameter
    addParameter( homogeneousMean );
}


void PhyloOrnsteinUhlenbeckProcess::setMean( const TypedDagNode< RbVector< double > >* x)
{
    // remove the old parameter
    removeParameter( heterogeneousMean );
    removeParameter( homogeneousMean );
    
    // set our local pointer copy
    heterogeneousMean = x;
    
    // add the new parameter
    addParameter( heterogeneousMean );
}


void PhyloOrnsteinUhlenbeckProcess::setPhi( const TypedDagNode< double >* x)
{
    // remove the old parameter
    removeParameter( heterogeneousPhi );
    removeParameter( homogeneousPhi );
    
    // set our local pointer copy
    homogeneousPhi = x;
    
    // add the new parameter
    addParameter( homogeneousPhi );
}


void PhyloOrnsteinUhlenbeckProcess::setPhi( const TypedDagNode< RbVector< double > >* x)
{
    // remove the old parameter
    removeParameter( heterogeneousPhi );
    removeParameter( homogeneousPhi );
    
    // set our local pointer copy
    heterogeneousPhi = x;
    
    // add the new parameter
    addParameter( heterogeneousPhi );
}


void PhyloOrnsteinUhlenbeckProcess::setSigma( const TypedDagNode< double >* x)
{
    // remove the old parameter
    removeParameter( heterogeneousSigma );
    removeParameter( homogeneousSigma );
    
    // set our local pointer copy
    homogeneousSigma = x;
    
    // add the new parameter
    addParameter( homogeneousSigma );
}
    
    
    
    
void PhyloOrnsteinUhlenbeckProcess::setSigma( const TypedDagNode< RbVector< double > >* x)
{
    // remove the old parameter
    removeParameter( heterogeneousSigma );
    removeParameter( homogeneousSigma );
        
    // set our local pointer copy
    heterogeneousSigma = x;
        
    // add the new parameter
    addParameter( heterogeneousSigma );
}


void PhyloOrnsteinUhlenbeckProcess::setRoot( const TypedDagNode< double >* x)
{
    // remove the old parameter
    removeParameter( root );
    
    // set our local pointer copy
    root = x;
    
    // add the new parameter
    addParameter( root );
}





