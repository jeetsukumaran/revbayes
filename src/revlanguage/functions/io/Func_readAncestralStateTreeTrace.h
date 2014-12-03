/**
 * @file
 * This file contains the declaration of Func_readAncestralStateTreeTrace.
 * This class is differentiated from Func_readTreeTrace only in the sense that
 * the trees put into the TreeTrace have not been rerooted; the nodes keep the
 * indexes already assigned to them, which is necessary for ancestral state
 * reconstruction.
 *
 * @brief Declaration of Func_readAncestralStateTreeTrace
 *
 * (c) Copyright 2014- under GPL version 3
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 */


#ifndef Func_readAncestralStateTreeTrace_H
#define Func_readAncestralStateTreeTrace_H

#include "RlBranchLengthTree.h"
#include "RlFunction.h"
#include "RbFileManager.h"
#include "RlTimeTree.h"
#include "RlTreeTrace.h"

#include <map>
#include <string>
#include <vector>


namespace RevLanguage {
    
    class Func_readAncestralStateTreeTrace :  public Function {
        
    public:
        // Basic utility functions
        Func_readAncestralStateTreeTrace*             clone(void) const;                                                                  //!< Clone the object
        static const std::string&       getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&          getClassTypeSpec(void);                                                             //!< Get class type spec
        const TypeSpec&                 getTypeSpec(void) const;                                                            //!< Get language type of the object
        
        // Regular functions
        RevPtr<Variable>                execute(void);                                                                      //!< Execute function
        const ArgumentRules&            getArgumentRules(void) const;                                                       //!< Get argument rules
        const TypeSpec&                 getReturnType(void) const;                                                          //!< Get type of return value
        
    private:
		
        void                            formatError(RevBayesCore::RbFileManager& fm, std::string& errorStr);                //!< Format the error string when (mis)reading files
        TreeTrace<BranchLengthTree>*    readBranchLengthTrees(const std::vector<std::string> &fns, const std::string &d);
        TreeTrace<TimeTree>*            readTimeTrees(const std::vector<std::string> &fns, const std::string &d);
    };
    
}

#endif

