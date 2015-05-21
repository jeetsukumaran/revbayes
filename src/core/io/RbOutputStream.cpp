#ifdef AP_MPI
#include <mpi.h>
#endif

#include "RbOutputStream.h"

using namespace RevBayesCore;

RbOutputStream::RbOutputStream( void )
{
    
}

//template<class T>
//RbOutputStream& RbOutputStream::operator<<(const T& x)
//{
//    
//    size_t pid = 0;
//#ifdef AP_MPI
//    pid = MPI::COMM_WORLD.Get_rank();
//#endif
//    if ( pid == 0 )
//    {
//        std::cout << x;
//    }
//    
//    return *this;
//}
