//
//  AbstractCharacterHistoryCtmc.cpp
//  rb_mlandis
//
//  Created by Michael Landis on 8/6/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#include "AbstractCharacterHistoryCtmc.h"
#include "DistributionExponential.h"
#include "ExponentialDistribution.h"
#include "TypedDistribution.h"
#include "RateMatrix.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"

using namespace RevBayesCore;

/*
AbstractCharacterHistoryCtmc::AbstractCharacterHistoryCtmc(BranchHistory* bh, TypedDagNode<RateMatrix> *rateMtx, size_t ns, size_t nc) : TypedDistribution<BranchHistory>(bh), rateMatrix(rateMtx), numStates(ns), numCharacters(nc)
{
    addParameter(rateMatrix);
}
*/

AbstractCharacterHistoryCtmc::AbstractCharacterHistoryCtmc(TypedDagNode<RateMatrix> *rateMtx, std::vector<const TypedDagNode<double>* > r, const TypedDagNode<TimeTree>* t, const TypedDagNode<double>* br, size_t nc, size_t ns, size_t idx) : TypedDistribution<BranchHistory>(new BranchHistory(nc, ns)), rateMatrix(rateMtx), numStates(ns), numCharacters(nc), rates(r), branchRate(br), tree(t), index(idx)
{
    for (size_t i = 0; i < r.size(); i++)
        addParameter(rates[i]);
    
    addParameter(rateMatrix);
    addParameter(tree);
    addParameter(branchRate);
    
    //redrawValue();
}

AbstractCharacterHistoryCtmc::AbstractCharacterHistoryCtmc(const AbstractCharacterHistoryCtmc& m) : TypedDistribution<BranchHistory>(m)
{
    rateMatrix = m.rateMatrix;
    rates = m.rates;
    tree = m.tree;
    branchRate = m.branchRate;
    numStates = m.numStates;
    numCharacters = m.numCharacters;
    index = m.index;
}

double AbstractCharacterHistoryCtmc::computeLnProbability(void)
{
    return 0.0;
}

void AbstractCharacterHistoryCtmc::simulatePath(void)
{
    value->clearEvents();
    
    // get transition rates
    std::vector<std::vector<double> > r(numStates);
    for (size_t i = 0; i < numStates; i++)
    {
        r[i].resize(numStates);
        for (size_t j = 0; j < numStates; j++)
        {
            if (i != j)
            {
                size_t k = i * (numStates-1) + j;
                r[i][i] -= rates[k]->getValue();
                r[i][j] = rates[k]->getValue();
            }
        }
    }

    // reject sample path history
    std::vector<CharacterEvent*> parentVector = value->getParentCharacters();
    std::vector<CharacterEvent*> childVector = value->getChildCharacters();
    std::set<CharacterEvent*,CharacterEventCompare> history;
    
    double br = branchRate->getValue();
    double bt = tree->getValue().getBranchLength(index);
    if (bt == 0.0)
        bt = 200;
    bt = 1.0;
    
    for (size_t i = 0; i < numCharacters; i++)
    {
        std::set<CharacterEvent*> tmpHistory;

        unsigned int currState = parentVector[i]->getState();
        unsigned int endState = childVector[i]->getState();
        do
        {
            tmpHistory.clear();
            double t = 0.0;
            do
            {
                double r_sum = -r[currState][currState];
                t += RbStatistics::Exponential::rv( r_sum * bt * br, *GLOBAL_RNG );
        
                if (t < 1.0)
                {
                    double u = GLOBAL_RNG->uniform01() * r_sum;
                    for (unsigned int j = 0; j < numStates; j++)
                    {
                        //std::cout << "sample " << i << " " << currState << " -> "<<  j << " " << u << " " << r_sum << "\n";
                        if (j != currState)
                        {
                            u -= r[currState][j];
                            if (u <= 0.0)
                            {
                                currState = j;
                                break;
                            }
                        }
                    }
                    CharacterEvent* evt = new CharacterEvent(i,currState,t);
                    tmpHistory.insert(evt);
                }
                //std::cout << i << " " << t << " " << currState << "\n";
            }
            while(t < 1.0);
            //std::cout << "--------\n";
        }
        while (currState != endState);
        //std::cout << "^ accepted\n--------\n";
        
        for (std::set<CharacterEvent*>::iterator it = tmpHistory.begin(); it != tmpHistory.end(); it++)
            history.insert(*it);
    }
    
    value->setHistory(history);
}

void AbstractCharacterHistoryCtmc::simulateChildCharacterState(void)
{
    value->setChildCharacters(simulateCharacterState(1.0));
}

void AbstractCharacterHistoryCtmc::simulateParentCharacterState(void)
{
    value->setParentCharacters(simulateCharacterState(0.0));
}

std::set<CharacterEvent*> AbstractCharacterHistoryCtmc::simulateCharacterState(double t)
{
    std::set<CharacterEvent*> characterState;

    for (size_t i = 0; i < numCharacters; i++)
    {
        unsigned int s = (unsigned int)(GLOBAL_RNG->uniform01() * numStates);
        //std::cout << s << " ";
        characterState.insert(new CharacterEvent(i,s,t));
    }
    //std::cout << "\n";
    
    return characterState;
}


void AbstractCharacterHistoryCtmc::redrawValue(void)
{
    
    //std::cout << "ns " << numStates << "   nc " << numCharacters << "\n";

    if (value->getRedrawParentCharacters())
    {
        //std::cout << index << " !! redraw parent\n";
        simulateParentCharacterState();
        value->setRedrawParentCharacters(false);
    }

    
    if (value->getRedrawChildCharacters())
    {
        //std::cout << index << " redraw child\n";
        simulateChildCharacterState();
        value->setRedrawChildCharacters(false);
    }
    

    if (value->getRedrawHistory())
    {
        //std::cout << index << " redraw path\n";
        simulatePath();
        value->setRedrawHistory(false);
    }
    
    //value->print();
}

void AbstractCharacterHistoryCtmc::swapParameter(const DagNode *oldP, const DagNode *newP)
{
    for (size_t i = 0; i < rates.size(); i++)
    {
        if (rates[i] == oldP)
        {
            rates[i] = static_cast<const TypedDagNode<double>* >(newP);
        }
    }
    if (tree == oldP)
    {
        tree = static_cast<const TypedDagNode<TimeTree>* >(newP);
    }
    else if (branchRate == oldP)
    {
        branchRate = static_cast<const TypedDagNode<double>* >(oldP);
    }
}

/*
void AbstractCharacterHistoryCtmc::setValue(BranchHistory* v)
{
    
}
*/


void AbstractCharacterHistoryCtmc::setValue(const std::multiset<CharacterEvent*,CharacterEventCompare>& updateSet, const std::set<CharacterEvent*>& parentSet, const std::set<CharacterEvent*>& childSet, const std::set<size_t>& indexSet )
{
    value->updateHistory(updateSet,parentSet,childSet,indexSet);
}


/*
 
 // silly distribution that samples according to rates away from some state
 std::vector<double> p;
 double p_sum = 0.0;
 for (size_t i = 0; i < numStates; i++)
 {
 double p_i = 0.0;
 for (size_t j = 0; j < numStates-1; j++)
 {
 if (i != j)
 {
 p_i += rates[i * numStates + j]->getValue();
 }
 }
 p.push_back(1.0 / p_i);
 p_sum += p[i];
 }
 
 for (size_t i = 0; i < p.size(); i++)
 p[i] /= p_sum;
 
 
 for (size_t i = 0; i < numCharacters; i++)
 {
 unsigned int s = 0;
 double u = GLOBAL_RNG->uniform01();
 for (size_t j = 0; j < p.size(); j++)
 {
 u -= p[j];
 if (u < 0.0)
 {
 s = 0;
 break;
 }
 }
 
 characterState.insert(new CharacterEvent(i,s,t));
 }
 */
