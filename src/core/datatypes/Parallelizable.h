#ifndef Parallelizable_H
#define Parallelizable_H

#include <string>

namespace RevBayesCore {
    
    
    /**
     * Interace for Parallelizable classes.
     *
     * The Parallelizable interface provides a mechanism for code parallelization.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-05-20, version 1.0
     */
    class Parallelizable {
        
    public:
        virtual                         ~Parallelizable(void) {}
        
        void                            setNumberOfProcesses(size_t n, size_t offset=0);                        //!< Set the number of processes.

    protected:
        
        // hidden constrocturs
        Parallelizable();                                                                                       //!< Constructor
        Parallelizable(const Parallelizable &n);                                                                //!< Copy Constructor
        
        Parallelizable&                 operator=(const Parallelizable &p);                                     //!< Overloaded assignment operator

        
        // protected methods that derived classes can overwrite
        virtual void                    setNumberOfProcessesSpecialized(size_t n, size_t offset=0);             //!< Set the number of processes in a specialized way for derived classes.
        

        // protected members available for derived classes
        size_t                          activePID;
        size_t                          numProcesses;
        size_t                          pid;
        bool                            processActive;

    };
    
}

#endif

