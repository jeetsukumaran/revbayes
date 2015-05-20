#include "Parallelizable.h"

using namespace RevBayesCore;


/**
 * Default constructor.
 * We initialize the variables used for code parallelization,
 * such as the number of processes, the process ID, etc.
 */
Parallelizable::Parallelizable() :
    activePID( 0 ),
    numProcesses( 1 ),
    pid( 0 ),
    processActive( true )
{
    
    
#ifdef RB_MPI
    //    numProcesses = MPI::COMM_WORLD.Get_size();
    pid = MPI::COMM_WORLD.Get_rank();
    processActive = (pid == activePID);
#endif
    
    
}
 

/**
 * Standard copy constructor.
 */
Parallelizable::Parallelizable( const Parallelizable &p ) :
    activePID( p.activePID ),
    numProcesses( p.numProcesses ),
    pid( p.pid ),
    processActive( p.processActive )
{
    
}


/**
 * Standard assignment operator.
 */
Parallelizable& Parallelizable::operator=(const Parallelizable &p)
{
    
    // check for self-assignment
    if ( this != &p )
    {
        activePID       = p.activePID;
        numProcesses    = p.numProcesses;
        pid             = p.pid;
        processActive   = p.processActive;
    }
    
    // return reference to myself
    return *this;
}


/**
 * Public method for setting the number of available processes and the process offset for this object.
 */
void Parallelizable::setNumberOfProcesses(size_t n, size_t offset)
{
    
    numProcesses    = n;
    activePID       = offset;
    processActive   = (pid == activePID);
    
    // delegate call for derived classes
}


/**
 * Dummy implementation which derived classes can overwrite.
*/
void Parallelizable::setNumberOfProcessesSpecialized(size_t n, size_t offset )
{
    
    // nothing done here.
}