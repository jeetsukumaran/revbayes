#ifndef Assignable_H
#define Assignable_H

#include <string>

namespace RevBayesCore {
    
    
    
    /**
     * Interace for assignable classes.
     *
     * The assignable interface provides a mechanism for safe inheritance and assignment
     * objects without knowing its actual derived typed.
     *
     * The only method this interface defines
     * is T& assign(const T&). The assign function is required in e.g. DAG nodes if the
     * type is abstract (e.g. RateMatrix, Function, ...) and for vectors. 
     * See also the Assign.h file.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-03-01, version 1.0
     */
    class Assignable {
        
    public:
        virtual                         ~Assignable(void) {}
        
        virtual Assignable&              assign( const Assignable &a ) = 0;         //!< Assign the object and prevent slicing
    };
    
}

#endif

