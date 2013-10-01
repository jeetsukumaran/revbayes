//
//  BranchHistory.cpp
//  rb_mlandis
//
//  Created by Michael Landis on 8/6/13.
//  Copyright (c) 2013 Michael Landis. All rights reserved.
//

#include "BranchHistory.h"
#include "CharacterEvent.h"
#include "CharacterEventCompare.h"

#include <iostream>
#include <iomanip>

using namespace RevBayesCore;


//BranchHistory::BranchHistory(void) : numCharacters(0), numStates(0), isTip(false), isRoot(false) { }

BranchHistory::BranchHistory(size_t nc, size_t ns) : numCharacters(nc), numStates(ns), redrawChildCharacters(true), redrawParentCharacters(true), redrawHistory(true) //, isTip(it), isRoot(ir)
{
    parentCharacters.resize(numCharacters);
    childCharacters.resize(numCharacters);
}


BranchHistory& BranchHistory::operator=(const BranchHistory &bh) {
    
    if (this != &bh) {
        
        this->numStates = bh.numStates;
        this->numCharacters = bh.numCharacters;
        this->parentCharacters = bh.parentCharacters;
        this->childCharacters = bh.childCharacters;
        this->history = bh.history;
        
        this->redrawChildCharacters = bh.redrawChildCharacters;
        this->redrawParentCharacters = bh.redrawParentCharacters;
        this->redrawHistory = bh.redrawHistory;
    }
    
    return *this;
}

BranchHistory* BranchHistory::clone(void) const
{
    return new BranchHistory(*this);
}

const std::vector<CharacterEvent*>& BranchHistory::getParentCharacters(void)
{
    return parentCharacters;
}

const std::vector<CharacterEvent*>& BranchHistory::getChildCharacters(void)
{
    return childCharacters;
}


const std::vector<CharacterEvent*>& BranchHistory::getParentCharacters(void) const
{
    return parentCharacters;
}

const std::vector<CharacterEvent*>& BranchHistory::getChildCharacters(void) const
{
    return childCharacters;
}

const std::multiset<CharacterEvent*,CharacterEventCompare>& BranchHistory::getHistory(void)
{
    return history;
}

void BranchHistory::clearEvents(void)
{
    history.clear();
}

void BranchHistory::clearEvents(const std::set<size_t>& indexSet)
{
    std::set<size_t>::iterator it_idx;
    std::multiset<CharacterEvent*,CharacterEventCompare>::iterator it_h, it_tmp;

    // store history elements that match indexSet into clearSet
    std::set<CharacterEvent*,CharacterEventCompare> clearSet;
    
    // for each event in history, delete if index matches indexSet
    for (it_h = history.begin(); it_h != history.end(); )
    {
        if ( indexSet.find( (*it_h)->getIndex() ) != indexSet.end() )
        {
            it_tmp = it_h;
            ++it_tmp;
            history.erase(it_h);
            it_h = it_tmp;
        }
        else
        {
            ++it_h;
        }
    }
}

void BranchHistory::addEvent(CharacterEvent* evt)
{
    history.insert(evt);
}

void BranchHistory::removeEvent(CharacterEvent* evt)
{
    history.erase(evt);
}

void BranchHistory::updateHistory(const std::multiset<CharacterEvent*,CharacterEventCompare>& updateSet, const std::set<CharacterEvent*>& parentSet, const std::set<CharacterEvent*>& childSet, const std::set<size_t>& indexSet)
{
    // erase events on branchHistory for indices in indexSet
    clearEvents(indexSet);
    
    // insert elements into history
    std::multiset<CharacterEvent*,CharacterEventCompare>::iterator it_h;
    for (it_h = updateSet.begin(); it_h != updateSet.end(); it_h++)
        history.insert(*it_h);
    
    // update events on terminal vectors
    std::set<CharacterEvent*>::iterator it_idx;
    for (it_idx = parentSet.begin(); it_idx != parentSet.end(); it_idx++)
        parentCharacters[ (*it_idx)->getIndex() ] = *it_idx;
    for (it_idx = childSet.begin(); it_idx != childSet.end(); it_idx++)
        childCharacters[ (*it_idx)->getIndex() ] = *it_idx;
    
}

void BranchHistory::setParentCharacters(const std::set<CharacterEvent*>& s)
{
    std::set<CharacterEvent*>::iterator it;
    for (it = s.begin(); it != s.end(); it++)
        parentCharacters[ (*it)->getIndex() ] = *it;
}

void BranchHistory::setChildCharacters(const std::set<CharacterEvent*>& s)
{
    std::set<CharacterEvent*>::iterator it;
    for (it = s.begin(); it != s.end(); it++)
        childCharacters[ (*it)->getIndex() ] = *it;
}


void BranchHistory::setHistory(const std::set<CharacterEvent*,CharacterEventCompare>& s)
{
    history.clear();
    for (std::set<CharacterEvent*,CharacterEventCompare>::iterator it = s.begin(); it != s.end(); it++)
        history.insert(*it);
}

void BranchHistory::setHistory(const std::multiset<CharacterEvent*,CharacterEventCompare>& s)
{
    history = s;
}


void BranchHistory::print(void)
{
    std::set<CharacterEvent*,CharacterEventCompare>::iterator it_h;
    std::vector<CharacterEvent*>::iterator it_v;
    std::vector<CharacterEvent*> tmp = parentCharacters;
    
    std::cout << "0.0      : ";
    for (it_v = parentCharacters.begin(); it_v != parentCharacters.end(); it_v++)
    {
        std::cout << std::setw(8) << (*it_v)->getState() << " ";
    }
    std::cout << "\n";
    
    for (it_h = history.begin(); it_h != history.end(); it_h++)
    {
        std::cout << (*it_h)->getTime() << " : ";
        tmp[ (*it_h)->getIndex() ] = *it_h;
        for (size_t i = 0; i < numCharacters; i++)
        {
            if (i != (*it_h)->getIndex())
                std::cout << std::setw(8) << " ";
            else
                std::cout << std::setw(8) << (*it_h)->getState();
            std::cout << " ";
        }
        std::cout << "\n";
    }
    std::cout << "1.0      : ";
    for (it_v = childCharacters.begin(); it_v != childCharacters.end(); it_v++)
    {
        std::cout << std::setw(8) << (*it_v)->getState() << " ";
    }
    std::cout << "\n";
    ;
}

const std::set<size_t>& BranchHistory::getDirtyCharacters(void)
{
    return dirtyCharacters;
}

void BranchHistory::setDirtyCharacters(const std::set<size_t>& s)
{
    dirtyCharacters = s;
}

bool BranchHistory::getRedrawChildCharacters(void)
{
    return redrawChildCharacters;
}

bool BranchHistory::getRedrawParentCharacters(void)
{
    return redrawParentCharacters;
}

bool BranchHistory::getRedrawHistory(void)
{
    return redrawHistory;
}

void BranchHistory::setRedrawChildCharacters(bool tf)
{
    redrawChildCharacters = tf;
}

void BranchHistory::setRedrawParentCharacters(bool tf)
{
    redrawParentCharacters = tf;
}

void BranchHistory::setRedrawHistory(bool tf)
{
    redrawHistory = tf;
}

std::ostream& RevBayesCore::operator<<(std::ostream& o, const BranchHistory& x) {
    
    o << " ";
    
    return o;
}