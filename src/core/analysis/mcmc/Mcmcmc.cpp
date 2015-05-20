#include "Mcmcmc.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RlUserInterface.h"
#include "RbConstants.h"
#include "RbException.h"

#include <iostream>
#include <vector>
#include <cmath>

#ifdef RB_MPI
#include <mpi.h>
#endif

using namespace RevBayesCore;

Mcmcmc::Mcmcmc(const Model& m, const RbVector<Move> &mv, const RbVector<Monitor> &mn, std::string sT, size_t nc, size_t si, double dt) : MonteCarloSampler( ),
    numChains(nc),
    scheduleType(sT),
    currentGeneration(0),
    swapInterval(si),
    activeChainIndex( 0 ),
    delta( dt ),
    generation( 0 ),
    numAttemptedSwaps( 0 ),
    numAcceptedSwaps( 0 )
{
    
//    // only use a many processes as we have chains
//    if (numChains < numProcesses)
//    {
//        numProcesses = numChains;
//    }
    
    // initialize container sizes
    chains.resize(numChains);
    chainValues.resize(numChains, 0.0);
    chainHeats.resize(numChains, 0.0);
    
    // assign chains to processors, instantiate Mcmc objects
    baseChain = new Mcmc(m, mv, mn);
    // assign chains to processors, instantiate Mcmc objects
    for (size_t i = 0; i < numChains; i++)
    {
        size_t process_index_start = size_t(floor( (double(i)   / numChains ) * numProcesses) );
        size_t process_index_end   = size_t(floor( (double(i+1) / numChains ) * numProcesses) );
        size_t num_processes       = process_index_end - process_index_start;
        
        // all chains know heat-order and chain-processor schedules
        heatRanks.push_back(i);
        
        // add chain to pid's chain vector (smaller memory footprint)
        if ( pid >= process_index_start && pid < process_index_end)
        {
            // get chain heat
            double b = computeBeta(delta, i);
            
            // create chains
            // @Michael: Why do we need to create new objects??? (Sebastian)
            Mcmc* oneChain = new Mcmc( *baseChain );
            oneChain->setScheduleType( scheduleType );
            oneChain->setChainActive( i == 0 );
            oneChain->setChainPosteriorHeat( b );
            oneChain->setChainIndex( i );
            oneChain->setNumberOfProcesses( num_processes, process_index_start);
            chains[i] = oneChain;
        }
        else
        {
            chains[i] = NULL;
        }
    }
    
}

Mcmcmc::Mcmcmc(const Mcmcmc &m) : MonteCarloSampler(m)
{
    
    delta               = m.delta;
    numChains           = m.numChains;
    numProcesses        = m.numProcesses;
    heatRanks           = m.heatRanks;
    swapInterval        = m.swapInterval;
    activeChainIndex    = m.activeChainIndex;
    scheduleType        = m.scheduleType;
    
    numAttemptedSwaps   = m.numAttemptedSwaps;
    numAcceptedSwaps    = m.numAcceptedSwaps;
    generation          = m.generation;
    
    chains.clear();
    chains.resize(numChains, NULL);
    for (size_t i = 0; i < m.chains.size(); i++)
    {
        if (m.chains[i] != NULL)
        {
            chains[i]   = m.chains[i]->clone();
        }
    }
    
    chainValues         = m.chainValues;
    chainHeats          = m.chainHeats;
    
    currentGeneration   = m.currentGeneration;
    baseChain           = m.baseChain->clone();
    
}

Mcmcmc::~Mcmcmc(void)
{
    for (size_t i = 0; i < chains.size(); i++)
    {
        if (chains[i] != NULL)
        {
            delete chains[i];
        }
    }
    chains.clear();
    delete baseChain;
}


Mcmcmc& Mcmcmc::operator=(const Mcmcmc &m)
{
    // check for self-assignment
    if ( this != &m )
    {
            
        delta               = m.delta;
        numChains           = m.numChains;
        numProcesses        = m.numProcesses;
        heatRanks           = m.heatRanks;
        swapInterval        = m.swapInterval;
        activeChainIndex    = m.activeChainIndex;
        scheduleType        = m.scheduleType;
            
        numAttemptedSwaps   = m.numAttemptedSwaps;
        numAcceptedSwaps    = m.numAcceptedSwaps;
        generation          = m.generation;
            
        chains.clear();
        chains.resize(numChains, NULL);
        for (size_t i = 0; i < m.chains.size(); i++)
        {
            if (m.chains[i] != NULL)
            {
                chains[i]   = m.chains[i]->clone();
            }
        }
            
        chainValues         = m.chainValues;
        chainHeats          = m.chainHeats;
        
        currentGeneration   = m.currentGeneration;
        baseChain           = m.baseChain->clone();
     
    }
    
    // return reference to myself
    return *this;
}


void Mcmcmc::initialize(void)
{
    
}


double Mcmcmc::computeBeta(double d, size_t idx)
{

    return 1.0 / (1.0+delta*idx);
}


Mcmcmc* Mcmcmc::clone(void) const
{
    return new Mcmcmc(*this);
}


double Mcmcmc::getModelLnProbability( void )
{
    synchronizeValues();
    
    for (size_t i=0; i<chains.size(); i++)
    {
        if ( chainHeats[i] == 1.0 )
        {
            return chainValues[i];
        }
    }
    
    return RbConstants::Double::neginf;
}

std::string Mcmcmc::getStrategyDescription( void ) const
{
    std::string description = "";
    std::stringstream stream;
    stream << "The MCMCMC simulator runs 1 cold chain and " << (numChains-1) << " heated chains.\n";
//    stream << chains[ chainsPerProcess[pid][0] ]->getStrategyDescription();
    description = stream.str();
    
    return description;
}


void Mcmcmc::initializeSampler( bool priorOnly )
{
    
    // initialize each chain
    for (size_t i = 0; i < chains.size(); i++)
    {
        if (chains[i] != NULL)
        {
            chains[i]->initializeSampler( priorOnly );
        }
    }
    
}


void Mcmcmc::monitor(unsigned long g)
{
    
    for (size_t i = 0; i < chains.size(); i++)
    {
        if (chains[i] != NULL && chains[i]->isChainActive() )
        {
            chains[i]->monitor(g);
        }
    }
    
}

void Mcmcmc::nextCycle(bool advanceCycle)
{
    
    // run each chain for this process
    for (size_t j = 0; j < chains.size(); j++)
    {
        
        if ( chains[j] != NULL )
        {
            // advance chain j by a single cycle
            chains[j]->nextCycle( advanceCycle );
        }
        
    } // loop over chains for this process
    
    if ( advanceCycle == true )
    {
        // advance gen counter
        ++currentGeneration;
    }
    
    if ( currentGeneration % swapInterval == 0 )
    {
        
        // perform chain swap
        for (size_t i = 0; i < numChains; i++)
        {
            swapChains();
        }
    }
    
}

void Mcmcmc::printOperatorSummary(void) const
{
    for (size_t i = 0; i < numChains; i++)
    {
        if (chains[i] != NULL)
        {
            chains[i]->printOperatorSummary();
        }
    }
}


void Mcmcmc::reset( void )
{
    
    // reset counters
    numAcceptedSwaps = 0;
    numAttemptedSwaps = 0;
    
//    /* Reset the monitors */
//    for (size_t i = 0; i < chainsPerProcess[pid].size(); i++)
//    {
//        RbVector<Monitor>& monitors = chains[ chainsPerProcess[pid][i] ]->getMonitors();
//        for (size_t i=0; i<monitors.size(); i++)
//        {
//            monitors[i].reset();
//        }
//    }
    
    for (size_t i = 0; i < chains.size(); i++)
    {
        if ( chains[i] != NULL )
        {
            chains[i]->reset();
        }
    }

    
}


/**
 * Set the heat of the likelihood of the current chain.
 * This heat is used in posterior posterior MCMC algorithms to
 * heat the likelihood
 * The heat is passed to the moves for the accept-reject mechanism.
 */
void Mcmcmc::setLikelihoodHeat(double h)
{
    
    for (size_t i = 0; i < chains.size(); i++)
    {
        if (chains[i] != NULL)
        {
            chains[i]->setLikelihoodHeat( h );
        }
    }
    
}

void Mcmcmc::setNumberOfProcessesSpecialized(size_t n, size_t offset)
{
    
    // @MJL: Note to self. The ctor assumes numProcesses==1, so all chains are assigned to that processor.
    // After cloning all chains across processors, you will then want to thin out the chains as needed.
    // This should behave much like the old Mcmcmc ctor code, except it cannot assume a fresh object state.
    
//    // only use a many processes as we have chains
//    if (numChains < numProcesses)
//    {
//        numProcesses = numChains;
//    }
    
    // initialize container sizes
    for (size_t i = 0; i < chains.size(); i++)
    {
        if (chains[i] != NULL)
            delete chains[i];
    }
    
    chains.clear();
    chainValues.clear();
    chainHeats.clear();
    
    chains.resize(numChains);
    chainValues.resize(numChains, 0.0);
    chainHeats.resize(numChains, 0.0);
    
    // assign chains to processors, instantiate Mcmc objects
    for (size_t i = 0; i < numChains; i++)
    {
        size_t process_index_start = size_t(floor( (double(i)   / numChains ) * numProcesses) );
        size_t process_index_end   = size_t(floor( (double(i+1) / numChains ) * numProcesses) );
        size_t num_processes       = process_index_end - process_index_start;
        
        // all chains know heat-order and chain-processor schedules
        heatRanks.push_back(i);
        
        // add chain to pid's chain vector (smaller memory footprint)
        if ( pid >= process_index_start && pid < process_index_end)
        {
            // get chain heat
            double b = computeBeta(delta, i);
            
            // create chains
            // @Michael: Why do we need to create new objects??? (Sebastian)
            Mcmc* oneChain = new Mcmc( *baseChain );
            oneChain->setScheduleType( scheduleType );
            oneChain->setChainActive( i == 0 );
            oneChain->setChainPosteriorHeat( b );
            oneChain->setChainIndex( i );
            oneChain->setNumberOfProcesses( num_processes, process_index_start);
            chains[i] = oneChain;
        }
        else
        {
            chains[i] = NULL;
        }
    }
    
}


void Mcmcmc::setReplicateIndex(size_t index)
{

    this->replicateIndex = index;
    for (size_t i = 0; i < chains.size(); i++)
    {
        if (chains[i] != NULL)
        {
            chains[i]->setReplicateIndex(index);
        }
    }
}


void Mcmcmc::setStoneIndex(size_t index)
{
    
    for (size_t i = 0; i < chains.size(); i++)
    {
        if (chains[i] != NULL)
        {
            chains[i]->setStoneIndex(index);
        }
    }
}



void Mcmcmc::startMonitors( void )
{
    // Monitor
    for (size_t i = 0; i < chains.size(); i++)
    {
        if (chains[i] != NULL)
        {
            chains[i]->startMonitors();
        }
    }
    
}


void Mcmcmc::startMonitors(size_t numCycles)
{
    
    // Monitor
    for (size_t i = 0; i < chains.size(); i++)
    {
        
        if (chains[i] != NULL)
        {
            chains[i]->startMonitors( numCycles );
            
        
            // monitor chain activeIndex only
            if (chains[i]->isChainActive() )
            {
                chains[i]->monitor(0);
            }
        }
        
    }
    
}


void Mcmcmc::synchronizeValues(void)
{
    
    // synchronize chain values
    double results[numChains];
    for (size_t j = 0; j < numChains; j++)
    {
        results[j] = 0.0;
    }
    for (size_t j = 0; j < chains.size(); j++)
    {
        if (chains[j] != NULL)
        {
            results[j] = chains[j]->getModelLnProbability();
        }
    }

#ifdef RB_MPI
    if (activePID != pid)
    {

        MPI::COMM_WORLD.Send(&results, int(numChains), MPI::DOUBLE, (int)activePID, 0);

    }

#endif
    
    if (activePID == pid)
    {
#ifdef RB_MPI
        for (size_t i = 0; i < numProcesses; i++)
        {
            // ignore self
            if (pid == i)
            {
                continue;
            }
            
            double tmp_results[numChains];
            for (size_t j = 0; j < numChains; j++)
            {
                tmp_results[j] = 0.0;
            }
            MPI::COMM_WORLD.Recv(&tmp_results, int(numChains), MPI::DOUBLE, (int)i, 0);

            for (size_t j = 0; j < chains.size(); j++)
            {
                if (chains[j] != NULL)
                {
                    results[j] = tmp_results[j];
                }
            }
        }
#endif
        for (size_t i = 0; i < chainValues.size(); i++)
        {
            chainValues[i] = results[i];
        }
    }
    
}

void Mcmcmc::synchronizeHeats(void)
{
    
    // synchronize heat values
    double heats[numChains];
    for (size_t j = 0; j < numChains; j++)
    {
        heats[j] = 0.0;
    }
    for (size_t j = 0; j < chains.size(); j++)
    {
        if (chains[j] != NULL)
        {
            heats[j] = chains[j]->getChainPosteriorHeat();
        }
    }
    
#ifdef RB_MPI
    // share the heats accross processes
    if (activePID != pid)
    {

        MPI::COMM_WORLD.Send(&heats, (int)numChains, MPI::DOUBLE, (int)activePID, 0);

    }
#endif
    
    if ( activePID == pid )
    {
#ifdef RB_MPI
        for (size_t i = offset; i < numProcesses+offset; ++i)
        {
            if (pid == i)
                continue;
            
            double tmp_heats[numChains];
            for (size_t j = 0; j < numChains; j++)
            {
                tmp_heats[j] = 0.0;
            }
            
            MPI::COMM_WORLD.Recv(&tmp_heats, (int)numChains, MPI::DOUBLE, (int)i, 0);
            
            for (size_t j = 0; j < chains.size(); j++)
            {
                size_t k = chainsPerProcess[i][j];
                
                heats[k] = tmp_heats[k];
            }
        }
#endif
        for (size_t i = 0; i < chainValues.size(); i++)
        {
            chainHeats[i] = heats[i];
        }
    }
    
}


// MJL: allow swapChains to take a swap function -- e.g. pairwise swap for 1..n-1
void Mcmcmc::swapChains(void)
{
    
    size_t numChains = chains.size();
    
    // exit if there is only one chain
    if (numChains < 2)
    {
        return;
    }
    
    // send all chain values to pid 0
    synchronizeValues();
    
    // send all chain heats to pid 0
    synchronizeHeats();

    // swap chains
    swapNeighborChains();
    swapRandomChains();

}

//
//
//void Mcmcmc::swapNeighborChains(void)
//{
//    
//    size_t numAccepted = 0;
//    double lnProposalRatio = 0.0;
//    
//    //for (size_t i = 1; i < numChains; i++)
//    for (size_t i = numChains-1; i > 0; i--)
//    {
//        // swap?
//        bool accept = false;
//        // swap adjacent chains
//        size_t j = 0;
//        size_t k = 0;
//        
//        if (processActive == true)
//        {
//            ++numAttemptedSwaps;
//            
//            j = heatRanks[i-1];
//            k = heatRanks[i];
//            
//            // compute exchange ratio
//            double bj = chainHeats[j];
//            double bk = chainHeats[k];
//            double lnPj = chainValues[j];
//            double lnPk = chainValues[k];
//            double lnR = bj * (lnPk - lnPj) + bk * (lnPj - lnPk) + lnProposalRatio;
//            
//            // determine whether we accept or reject the chain swap
//            double u = GLOBAL_RNG->uniform01();
//            if (lnR >= 0)
//            {
//                accept = true;
//            }
//            else if (lnR < -100)
//            {
//                accept = false;
//            }
//            else if (u < exp(lnR))
//            {
//                accept = true;
//            }
//            else
//            {
//                accept = false;
//            }
//            
//            if (accept == true)
//            {
//                numAccepted++;
//            }
//            
//            // on accept, swap beta values and active chains
//            if (accept)
//            {
//                
//                //size_t tmpIdx = j;
//                heatRanks[i-1] = k;
//                heatRanks[i] = j;
//                
//                // swap active chain
//                if (activeChainIndex == j)
//                {
//                    activeChainIndex = k;
//                }
//                else if (activeChainIndex == k)
//                {
//                    activeChainIndex = j;
//                }
//                
//            }
//        }
//        
//        if (accept)
//        {
//            updateChainState(j);
//            updateChainState(k);
//            
//            ++numAcceptedSwaps;
//        }
//    }
//    
//}



void Mcmcmc::swapNeighborChains(void)
{
    
    double lnProposalRatio = 0.0;
    
    // randomly pick the indices of two chains
    int j = 0;
    int k = 0;
    
    // swap?
    bool accept = false;
    if (numChains < 2) return;
    
    if ( pid == activePID )
    {
        j = int(GLOBAL_RNG->uniform01() * (numChains-1));
        k = j + 1;
//        if (numChains > 1)
//        {
//            do {
//                k = int(GLOBAL_RNG->uniform01() * numChains);
//            }
//            while(j == k);
//        }
        
        ++numAttemptedSwaps;
        
        // compute exchange ratio
        double bj = chainHeats[j];
        double bk = chainHeats[k];
        double lnPj = chainValues[j];
        double lnPk = chainValues[k];
        double lnR = bj * (lnPk - lnPj) + bk * (lnPj - lnPk) + lnProposalRatio;
        
        // determine whether we accept or reject the chain swap
        double u = GLOBAL_RNG->uniform01();
        if (lnR >= 0)
        {
            accept = true;
        }
        else if (lnR < -100)
        {
            accept = false;
        }
        else if (u < exp(lnR))
        {
            accept = true;
        }
        else
        {
            accept = false;
        }
        
        
        // on accept, swap beta values and active chains
        if (accept == true )
        {
            
            // swap active chain
            if (activeChainIndex == j)
            {
                activeChainIndex = k;
            }
            else if (activeChainIndex == k)
            {
                activeChainIndex = j;
            }
            
            chainHeats[j] = bk;
            chainHeats[k] = bj;
            size_t tmp = heatRanks[j];
            heatRanks[j] = heatRanks[k];
            heatRanks[k] = tmp;
            
            ++numAcceptedSwaps;
        }
        
        
    }
    
#ifdef RB_MPI
    
    MPI::COMM_WORLD.Bcast(&j, 1, MPI_INT, (int)activePID);
    MPI::COMM_WORLD.Bcast(&k, 1, MPI_INT, (int)activePID);

#endif
    
    
    // update the chains accross processes
    // this is necessary because only process 0 does the swap
    // all the other processes need to be told that there was a swap
    updateChainState(j);
    updateChainState(k);
    
}



void Mcmcmc::swapRandomChains(void)
{
    
    double lnProposalRatio = 0.0;
    
    // randomly pick the indices of two chains
    int j = 0;
    int k = 0;
    
    // swap?
    bool accept = false;
    
    if ( pid == activePID )
    {
        j = int(GLOBAL_RNG->uniform01() * numChains);
        if (numChains > 1)
        {
            do {
                k = int(GLOBAL_RNG->uniform01() * numChains);
            }
            while(j == k);
        }
        
        ++numAttemptedSwaps;
            
        // compute exchange ratio
        double bj = chainHeats[j];
        double bk = chainHeats[k];
        double lnPj = chainValues[j];
        double lnPk = chainValues[k];
        double lnR = bj * (lnPk - lnPj) + bk * (lnPj - lnPk) + lnProposalRatio;
            
        // determine whether we accept or reject the chain swap
        double u = GLOBAL_RNG->uniform01();
        if (lnR >= 0)
        {
            accept = true;
        }
        else if (lnR < -100)
        {
            accept = false;
        }
        else if (u < exp(lnR))
        {
            accept = true;
        }
        else
        {
            accept = false;
        }
        
        
        // on accept, swap beta values and active chains
        if (accept == true )
        {
            
            // swap active chain
            if (activeChainIndex == j)
            {
                activeChainIndex = k;
            }
            else if (activeChainIndex == k)
            {
                activeChainIndex = j;
            }
            
            chainHeats[j] = bk;
            chainHeats[k] = bj;
            size_t tmp = heatRanks[j];
            heatRanks[j] = heatRanks[k];
            heatRanks[k] = tmp;
            
            ++numAcceptedSwaps;
        }
        
        
    }

#ifdef RB_MPI
    
    MPI::COMM_WORLD.Bcast(&j, 1, MPI_INT, (int)activePID);
    MPI::COMM_WORLD.Bcast(&k, 1, MPI_INT, (int)activePID);

#endif
    
    
    // update the chains accross processes
    // this is necessary because only process 0 does the swap
    // all the other processes need to be told that there was a swap
    updateChainState(j);
    updateChainState(k);
    
}


void Mcmcmc::tune( void )
{
    
    double rate = numAcceptedSwaps / double(numAttemptedSwaps);
    
    if ( rate > 0.44 )
    {
        delta *= (1.0 + ((rate-0.44)/0.56) );
    }
    else
    {
        delta /= (2.0 - rate/0.44 );
    }
    
}


void Mcmcmc::updateChainState(size_t j)
{
    
#ifdef RB_MPI
    MPI::COMM_WORLD.Barrier();
    // update heat
    if (pid == activePID && pid == processPerChain[j])
    {
        ; // do nothing
    }
    else if (pid == activePID)
    {
        MPI::COMM_WORLD.Send(&chainHeats[j], 1, MPI::DOUBLE, (int)processPerChain[j], 0);
    }
    else if (pid == processPerChain[j])
    {
        MPI::COMM_WORLD.Recv(&chainHeats[j], 1, MPI::DOUBLE, (int)activePID, 0);
    }
#endif
    
    if (chains[j] != NULL)
    {
        chains[j]->setChainPosteriorHeat(chainHeats[j]);
    }
    // update active state
    bool tf = activeChainIndex == j;
    
#ifdef RB_MPI
    if (pid == activePID && pid == processPerChain[j])
    {
        ; // do nothing
    }
    else if (pid == activePID)
    {
        MPI::COMM_WORLD.Send(&tf, 1, MPI::BOOL, (int)processPerChain[j], 0);
    }
    else if (pid == processPerChain[j])
    {
        MPI::COMM_WORLD.Recv(&tf, 1, MPI::BOOL, (int)activePID, 0);
    }
#endif
    
    if (chains[j] != NULL)
    {
        chains[j]->setChainActive( chainHeats[j] == 1.0 );
    }
    
}

