#ifndef Cloneable_H
#define Cloneable_H

#include <string>

namespace RevBayesCore {
    
    
    /**
     * Interace for cloneable classes.
     *
     * The cloneable interface provides a mechanism for safe inheritance and copying (cloning)
     * objects without knowing its actual derived typed.
     *
     * The only method this interface defines
     * is T* clone(void) const. The clone function is required in e.g. DAG nodes if the
     * type is abstract (e.g. RateMatrix, Function, ...)
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-03-01, version 1.0
     */
    class Cloneable {
        
    public:
        virtual                         ~Cloneable(void) {}
        
        virtual Cloneable*              clone( void ) const = 0;                                    //!< Create a clone/copy of the object
    };
    
}

#endif

