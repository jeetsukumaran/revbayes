//
//  RlRateGenerator.h
//  revbayes-proj
//
//  Created by Michael Landis on 3/17/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#ifndef __revbayes_proj__RlRateGenerator__
#define __revbayes_proj__RlRateGenerator__

#include "ModelObject.h"
#include "RateGenerator.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class RateGenerator : public ModelObject<RevBayesCore::RateGenerator> {
        
    public:
        
        RateGenerator(void);                                                                                                           //!< Default constructor
        RateGenerator(const RevBayesCore::RateGenerator& m);                                                                              //!< Default constructor
        RateGenerator(RevBayesCore::RateGenerator *m);                                                                                    //!< Default constructor
        RateGenerator(RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator> *d);                                                                                                        //!< Default constructor
        
        // Basic utility functions
        RateGenerator*                      clone(void) const;                                                                      //!< Clone object
        static const std::string&           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                                                 //!< Get class type spec
        const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        
        // Member method functions
        virtual RevPtr<RevVariable>         executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Map member methods to internal functions
        
    };
    
}

#endif /* defined(__revbayes_proj__RlRateGenerator__) */
