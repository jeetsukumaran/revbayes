/**
 * @file
 * This file contains the declaration of ArgumentRule, which is
 * the base class for objects used to describe rules for
 * arguments passed to functions.
 *
 * @brief Declaration of ArgumentRule
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-11-20, version 1.0
 *
 * $Id$
 */

#ifndef ArgumentRule_H
#define ArgumentRule_H

#include "RbInternal.h"
#include "TypeSpec.h"
#include "Environment.h"
#include "VariableSlot.h"

#include <string>

class DAGNode;
class RbObject;
class VectorString;

const std::string ArgumentRule_name = "argument rule";

class ArgumentRule : public RbInternal {

    public:
        // Basic utility functions
        virtual ArgumentRule*       clone(void) const { return new ArgumentRule(*this); }                                               //!< Clone object
        virtual const VectorString& getClass(void) const;                                                                               //!< Get class vector
        virtual const TypeSpec&     getTypeSpec(void) const;                                                                            //!< Get language type of the object
        void                        printValue(std::ostream& o) const;                                                                  //!< Print value for user
        std::string                 richInfo(void) const;                                                                               //!< General info on object

        // ArgumentRule functions
        const std::string&          getArgumentLabel(void) const;                                                                       //!< Get label of argument
        const std::string&          getArgumentType(void) const;                                                                        //!< Get argument type
        const TypeSpec&             getArgumentTypeSpec(void) const;                                                                    //!< Get argument type spec
        RbPtr<const Variable>       getDefaultVariable(void) const;                                                                     //!< Get default argument
        RbPtr<Variable>             getDefaultVariable(void);                                                                           //!< Get default argument (non-const to return non-const variable)
        bool                        hasDefault(void) const;                                                                             //!< Has default?
        virtual bool                isArgumentValid(RbPtr<const DAGNode> var, bool& needsConversion) const;                             //!< Is var valid argument?

    protected:
                                    ArgumentRule(const std::string& argName, RbPtr<RbLanguageObject> defValue);                         //!< Constructor of rule from default value
                                    ArgumentRule(const std::string& argName, const TypeSpec& argTypeSp);                                //!< Constructor of rule without default value
                                    ArgumentRule(const std::string& argName, const TypeSpec& argTypeSp, RbPtr<RbLanguageObject> defValue); //!< Constructor of rule with default value
                                    ArgumentRule(const std::string& argName, const TypeSpec& argTypeSp, RbPtr<DAGNode> defVariable);    //!< Constructor of rule with default reference or default wrapped value

        std::string                 label;                                                                                              //!< Label of argument
        VariableSlot                argSlot;                                                                                            //!< Slot with typespec and possibly default value
        bool                        hasDefaultVal;                                                                                      //!< Has default (which can be NULL) ?
    
    private:
        static const TypeSpec       typeSpec;
};

#endif
