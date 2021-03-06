#ifndef Move_BetaSimplex_H
#define Move_BetaSimplex_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Rev Wrapper of a beta-simplex move on a single elements of simplex.
     *
     * This class is the RevLanguage wrapper of ElementScale.
     *
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @copyright GPL version 3
     * @since 2015-05-21, version 1.0
     */
    class Move_BetaSimplex : public Move {
        
    public:
        
        Move_BetaSimplex(void);                                                                               //!< Default constructor (0.0)
        
        // Basic utility functions
        virtual Move_BetaSimplex*                   clone(void) const;                                                      //!< Clone object
        void                                        constructInternalObject(void);                                          //!< We construct the a new internal SlidingMove.
        static const std::string&                   getClassType(void);                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                 //!< Get class type spec
        const MemberRules&                          getParameterRules(void) const;                                             //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                      //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);     //!< Set member variable
        
        RevPtr<const RevVariable>                   x;
        RevPtr<const RevVariable>                   alpha;
        RevPtr<const RevVariable>                   tune;
        
    };
    
}

#endif
