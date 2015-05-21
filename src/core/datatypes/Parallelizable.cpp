#include "Parallelizable.h"


#ifdef RB_MPI
#include <mpi.h>
#endif

using namespace RevBayesCore;


/**
 * Default constructor.
 * We initialize the variables used for code parallelization,
 * such as the number of processes, the process ID, etc.
 */
Parallelizable::Parallelizable() :
    offsetPID( 0 ),
    numProcesses( 1 ),
    pid( 0 )
{
    
    
#ifdef RB_MPI
    numProcesses = MPI::COMM_WORLD.Get_size();
    pid = MPI::COMM_WORLD.Get_rank();
#endif
    
    
}
 

/**
 * Standard copy constructor.
 */
Parallelizable::Parallelizable( const Parallelizable &p ) :
    offsetPID( p.offsetPID ),
    numProcesses( p.numProcesses ),
    pid( p.pid )
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
        offsetPID       = p.offsetPID;
        numProcesses    = p.numProcesses;
        pid             = p.pid;
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
    offsetPID       = offset;
    
    // delegate call for derived classes
    setNumberOfProcessesSpecialized(n, offset);
}


/**
 * Dummy implementation which derived classes can overwrite.
*/
void Parallelizable::setNumberOfProcessesSpecialized(size_t n, size_t offset )
{
    
    // nothing done here.
}