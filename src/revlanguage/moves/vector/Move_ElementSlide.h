#ifndef Move_ElementSlide_H
#define Move_ElementSlide_H

#include "SimpleMove.h"
#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Rev Wrapper of a scaling move on all elements of a real valued vector.
     *
     * This class is the RevLanguage wrapper of ElementSlide.
     *
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @copyright GPL version 3
     * @since 2015-05-21, version 1.0
     */
    class Move_ElementSlide : public Move {
        
    public:
        
        Move_ElementSlide(void);                                                                                                                   //!< Default constructor
        
        // Basic utility functions
        virtual Move_ElementSlide*                  clone(void) const;                                                                              //!< Clone the object
        void                                        constructInternalObject(void);                                                                  //!< We construct the a new internal move.
        static const std::string&                   getClassType(void);                                                                             //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                         //!< Get class type spec
        const MemberRules&                          getParameterRules(void) const;                                                                     //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                        //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                              //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);             //!< Set member variable
        
        RevPtr<const RevVariable>                   x;                                                                                              //!< The variable holding the real valued vector.
        RevPtr<const RevVariable>                   delta;                                                                                         //!< The variable for the tuning parameter.
        RevPtr<const RevVariable>                   tune;                                                                                           //!< The variable telling if to tune or not.
        
    };
    
}


#endif /* defined(__revbayes_proj__Move_ElementSlide__) */
