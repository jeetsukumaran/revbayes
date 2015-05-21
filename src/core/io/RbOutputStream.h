#ifndef RbOutputStream_H
#define RbOutputStream_H

#include <iostream>
#include <ostream>

#ifdef AP_MPI
#include <mpi.h>
#endif

namespace RevBayesCore {
    
    
    
    /**
     * General output stream for RevBayes.
     *
     * This output stream prints to std::cout. It is for convenience so that it can be used with MPI.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-05-20, version 1.0
     *
     */
    class RbOutputStream : public std::ostream {
    
    public:
        RbOutputStream();

        template<class T>
        RbOutputStream& operator<<(const T& x)
        {
            size_t pid = 0;
#ifdef AP_MPI
            pid = MPI::COMM_WORLD.Get_rank();
#endif
            if ( pid == 0 )
            {
                std::cout << x;
            }
            
            return *this;
        }
        
    
    };
    
}

#endif
