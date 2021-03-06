#ifndef Move_UpDownTreeScale_H
#define Move_UpDownTreeScale_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the up-down-scaling move.
     *
     * The RevLanguage wrapper of the up-down-scaling move simply
     * manages the interactions through the Rev with our core.
     * That is, the internal move object can be constructed and hooked up
     * in a DAG-nove (variable) that it works on.
     * See the UpDownTreeScaleMove.h for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-05-25, version 1.0
     *
     */
    class Move_UpDownTreeScale : public Move {
        
    public:
        
        Move_UpDownTreeScale(void);                                                                                                            //!< Default constructor
        
        // Basic utility functions
        virtual Move_UpDownTreeScale*                   clone(void) const;                                                                          //!< Clone object
        void                                        constructInternalObject(void);                                                              //!< We construct the a new internal SlidingMove.
        static const std::string&                   getClassType(void);                                                                         //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                     //!< Get class type spec
        const MemberRules&                          getParameterRules(void) const;                                                                 //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                    //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                          //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);          //!< Set member variable
        
        RevPtr<const RevVariable>                   upScalar;                                                                                          //!< The variable on which the move works
        RevPtr<const RevVariable>                   upVector;                                                                                          //!< The variable on which the move works
        RevPtr<const RevVariable>                   upTree;                                                                                          //!< The variable on which the move works
        RevPtr<const RevVariable>                   downScalar;                                                                                           //!< The variable on which the move works
        RevPtr<const RevVariable>                   downVector;                                                                                           //!< The variable on which the move works
        RevPtr<const RevVariable>                   downTree;                                                                                           //!< The variable on which the move works
        RevPtr<const RevVariable>                   lambda;                                                                                          //!< The tuning parameter
        RevPtr<const RevVariable>                   tune;                                                                                      //!< The variable on which the move works
        
    };
    
}

#endif
