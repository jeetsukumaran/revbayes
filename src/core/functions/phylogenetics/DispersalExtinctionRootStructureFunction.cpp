//
//  DispersalExtinctionRootStructureFunction.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 3/3/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#include "DispersalExtinctionRootStructureFunction.h"
#include "RbException.h"

#include <cmath>
#include <vector>
#include <map>

using namespace RevBayesCore;

DispersalExtinctionRootStructureFunction::DispersalExtinctionRootStructureFunction(TypedDagNode<RbVector<double> >* rf, TypedDagNode<RbVector<double> >* rs) : TypedFunction<RbVector<double> >( new RbVector<double>() ),
rootFrequencies( rf ),
rangeSize(rs)
{
    // add the lambda parameter as a parent
    addParameter( rootFrequencies );
    addParameter( rangeSize );

    numStates = rootFrequencies->getValue().size();
    numCharacters = 0;
    size_t n = numStates;
    
    do {
        n /= 2;
        numCharacters++;
    } while(n != 1);

    (*value).resize(numStates, 0.0);
    
    
    makeBits();
    makeIdxByRangeSize();
   
    update();
}


DispersalExtinctionRootStructureFunction::~DispersalExtinctionRootStructureFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



DispersalExtinctionRootStructureFunction* DispersalExtinctionRootStructureFunction::clone( void ) const
{
    return new DispersalExtinctionRootStructureFunction( *this );
}


void DispersalExtinctionRootStructureFunction::keep(DagNode *affecter)
{
    //delegate to base class
    TypedFunction< RbVector<double> >::keep( affecter );
    
}

void DispersalExtinctionRootStructureFunction::makeBits(void)
{
    bits = std::vector<std::vector<unsigned> >(numStates, std::vector<unsigned>(numCharacters, 0));
    for (size_t i = 1; i < numStates; i++)
    {
        size_t n = i;
        for (size_t j = 0; j < numCharacters; j++)
        {
            bits[i][j] = n % 2;
            n /= 2;
            if (n == 0)
                break;
        }
    }
    for (size_t i = 0; i < (size_t)numStates; i++)
    {
        inverseBits[ bits[i] ] = (unsigned)i;
    }
}

void DispersalExtinctionRootStructureFunction::makeIdxByRangeSize(void)
{
    idxByRangeSize.resize(numCharacters + 1);
    for (size_t i = 0; i < numStates; i++)
    {
        int n = 0;
        for (size_t j = 0; j < bits[i].size(); j++)
            n += bits[i][j];
        idxByRangeSize[n].push_back((unsigned)i);
    }
}

void DispersalExtinctionRootStructureFunction::reInitialized( void )
{
    ;//    *value = tau->getValue();
}


void DispersalExtinctionRootStructureFunction::restore(DagNode *restorer)
{
    //delegate to base class
    TypedFunction< RbVector<double> >::restore( restorer );
}


void DispersalExtinctionRootStructureFunction::touch(DagNode *toucher)
{
    
    //delegate to base class
    TypedFunction< RbVector<double> >::touch( toucher );
    
}


void DispersalExtinctionRootStructureFunction::update( void )
{
    const RbVector<double>& rf = rootFrequencies->getValue();
    size_t minRangeSize = 1;
    
    // get normalized range sizes, s.t. largest range size value = 1
    const std::vector<double>& rs = rangeSize->getValue();
    double largestRangeSize = rs[0];
    for (size_t i = minRangeSize; i <= numCharacters; i++)
    {
        if (rs[i] > largestRangeSize)
        {
            largestRangeSize = rs[i];
        }
    }
    std::vector<double> normalizedRangeSize = rs;
    for (size_t i = 0; i < rs.size(); i++)
    {
        normalizedRangeSize[i] = rs[i] / largestRangeSize;
    }
    
    
    std::vector<std::vector<unsigned> >::iterator it;
    std::vector<unsigned>::iterator jt;
    
    *value = RbVector<double>(numStates, 0.0);
    
    double sum = 0.0;
    for (size_t i = minRangeSize; i <= numCharacters; i++)
    {
        for (std::vector<unsigned>::iterator it = idxByRangeSize[i].begin(); it != idxByRangeSize[i].end(); it++)
        {
            (*value)[ *it ] = rf[ *it ] * normalizedRangeSize[i];
            sum += rf[ *it ];
        }
    }
    for (size_t i = minRangeSize; i <= numCharacters; i++)
    {
        for (std::vector<unsigned>::iterator it = idxByRangeSize[i].begin(); it != idxByRangeSize[i].end(); it++)
        {
            (*value)[ *it ] = (*value)[ *it ] / sum;
        }
    }
}



void DispersalExtinctionRootStructureFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == rootFrequencies)
    {
        rootFrequencies = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    else if (oldP == rangeSize)
    {
        rangeSize = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
}
