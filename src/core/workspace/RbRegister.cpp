/**
 * @file
 * This file contains the Workspace function that adds types and functions
 * to the global workspace, registering them with the interpreter/compiler
 * during the process.
 *
 * @brief Function registering language objects
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @extends Frame
 * @package parser
 * @version 1.0
 * @since version 1.0 2009-09-02
 *
 * $Id:$
 */

#include "ConstantNode.h"
#include "ConstructorFunction.h"
#include "Container.h"
#include "Distribution.h"
#include "DistributionFunction.h"
#include "DistributionInterval.h"
#include "FunctionTable.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbAbstract.h"
#include "RbException.h"
#include "RbFunction.h"
#include "RbObject.h"
#include "UserInterface.h"
#include "ValueContainer.h"
#include "VariableContainer.h"
#include "Workspace.h"


/* Primitive types (alphabetic order) */
#include "AminoAcidState.h"
#include "Boolean.h"
#include "CharacterObservationContinuous.h"
#include "DnaState.h"
#include "Integer.h"
#include "Natural.h"
#include "RateMatrix.h"
#include "RbString.h"
#include "Real.h"
#include "RealPos.h"
#include "RnaState.h"
#include "StandardState.h"

/* Container types (alphabetic order) */
#include "MatrixReal.h"
#include "VectorAminoAcidStates.h"
#include "VectorBoolean.h"
#include "VectorContinuousObservations.h"
#include "VectorDnaStates.h"
#include "VectorInteger.h"
#include "VectorNatural.h"
#include "VectorReal.h"
#include "VectorRealPos.h"
#include "VectorRnaStates.h"
#include "VectorStandardStates.h"
#include "VectorString.h"

/* MemberObject types without auto-generated constructors(alphabetic order) */
#include "CharacterMatrix.h"
#include "List.h"
#include "Model.h"
#include "Simplex.h"

/* MemberObject types with auto-generated constructors (alphabetic order) */
#include "Mcmc.h"
#include "Move_mmultinomial.h"
#include "Move_mscale.h"
#include "Move_msimplex.h"
#include "Move_mslide.h"

/* Distributions with distribution constructors and distribution functions (alphabetic order) */
#include "Dist_dirichlet.h"
#include "Dist_exp.h"
#include "Dist_multinomial.h"
#include "Dist_norm.h"
#include "Dist_unif.h"

/* Basic internal functions (alphabetic order) */
#include "Func__range.h"
       
/* Basic arithmetic/logic templated functions */
#include "Func__add.h"
#include "Func__and.h"
#include "Func__div.h"
#include "Func__eq.h"
#include "Func__ge.h"
#include "Func__gt.h"
#include "Func__le.h"
#include "Func__lt.h"
#include "Func__mul.h"
#include "Func__ne.h"
#include "Func__or.h"
#include "Func__range.h"
#include "Func__sub.h"

/* Regular functions */
#include "Func_clamp.h"
#include "Func_ls.h"
#include "Func_model.h"
#include "Func_normalize.h"
#include "Func_quit.h"
#include "Func_s.h"
#include "Func_str.h"
#include "Func_v.h"

/* Math functions */
#include "Func_exp.h"
#include "Func_ln.h"
#include "Func_log.h"
#include "Func_power.h"
#include "Func_sin.h"
#include "Func_sqrt.h"

/* Phylogeny functions */
#include "Func_gtr.h"
#include "Func_readAlignment.h"

/* Regular templated functions */
#include "Func_transpose.h"

#include <sstream>
#include <vector>


/** Initialize global workspace */
void Workspace::initializeGlobalWorkspace(void) {

    try {
        /* Add types: add a dummy variable which we use for type checking, conversion checking and other tasks. */

        /* Add special abstract types that do not correspond directly to classes */
        addType( new RbAbstract( VectorString(RbVoid_name) ) );
        addType( new RbAbstract( VectorString(Scalar_name) ) );

        /* Add normal abstract types */
        addType( new RbAbstract( VectorString(RbObject_name) ) );
        addType( new RbAbstract( VectorString(Vector_name) + RbString(Container_name) + RbString(RbComplex_name) + RbString(RbObject_name), 1, Scalar_name ) );
        addType( new RbAbstract( VectorString(Matrix_name) + RbString(Container_name) + RbString(RbComplex_name) + RbString(RbObject_name), 2, Scalar_name ) );
        addType( new RbAbstract( VectorString(MemberObject_name) + RbString(RbComplex_name) + RbString(RbObject_name) ) );

        /* Add primitive types (alphabetic order) */
        addType( new AminoAcidState()                 );
        addType( new Boolean()                        );
        addType( new CharacterObservationContinuous() );
        addType( new DnaState()                       );
        addType( new Integer()                        );
        addType( new Natural()                        );
        addType( new RateMatrix()                     );
        addType( new RbString()                       );
        addType( new Real()                           );
        addType( new RealPos()                        );
        addType( new RnaState()                       );
        addType( new StandardState()                  );

        /* Add container types (alphabetic order) */
        addType( new MatrixReal()                   );
        addType( new VectorAminoAcidStates()        );
        addType( new VectorBoolean()                );
        addType( new VectorContinuousObservations() );
        addType( new VectorDnaStates()              );
        addType( new VectorInteger()                );
        addType( new VectorNatural()                );
        addType( new VectorReal()                   );
        addType( new VectorRealPos()                );
        addType( new VectorRnaStates()              );
        addType( new VectorStandardStates()         );
        addType( new VectorString()                 );

        /* Add MemberObject types without auto-generated constructors (alphabetic order) */
        addType( new CharacterMatrix() );
        addType( new List()            );
        addType( new Model()           );
        addType( new Simplex()         );

        /* Add MemberObject types with auto-generated constructors (alphabetic order) */
        addTypeWithConstructor( "mcmc",         new Mcmc()              );
        addTypeWithConstructor( "mmultinomial", new Move_mmultinomial() );
        addTypeWithConstructor( "mslide",       new Move_mslide()       );
        addTypeWithConstructor( "mscale",       new Move_mscale()       );
        addTypeWithConstructor( "msimplex",     new Move_msimplex()     );

        /* Add Distribution types with auto-generated constructors and distribution functions (alphabetic order) */
        addDistribution( "dirichlet",   new Dist_dirichlet()   );
        addDistribution( "exp",         new Dist_exp()         );
        addDistribution( "multinomial", new Dist_multinomial() );
        addDistribution( "norm",        new Dist_norm()        );
        addDistribution( "unif",        new Dist_unif()        );

        /* Now we have added all primitive and complex data types and can start type checking */
        Workspace::globalWorkspace().typesInitialized = true;
        Workspace::userWorkspace().typesInitialized   = true;

        /* Add basic internal functions (alphabetic order) */
        addFunction( "_range",    new Func__range()       );
       
        /* Add basic arithmetic templated functions */
        addFunction( "_add",      new Func__add<            Integer,        Integer,    Integer >() );
        addFunction( "_add",      new Func__add<               Real,           Real,       Real >() );
        addFunction( "_add",      new Func__add<            Integer,           Real,       Real >() );
        addFunction( "_add",      new Func__add<               Real,        Integer,       Real >() );
        addFunction( "_add",      new Func__add<         MatrixReal,     MatrixReal, MatrixReal >() );
        addFunction( "_div",      new Func__div<            Integer,        Integer,       Real >() );
        addFunction( "_div",      new Func__div<               Real,           Real,       Real >() );
        addFunction( "_div",      new Func__div<            Integer,           Real,       Real >() );
        addFunction( "_div",      new Func__div<               Real,        Integer,       Real >() );
        addFunction( "_div",      new Func__div<         MatrixReal,     MatrixReal, MatrixReal >() );
        addFunction( "_mul",      new Func__mul<            Integer,        Integer,    Integer >() );
        addFunction( "_mul",      new Func__mul<               Real,           Real,       Real >() );
        addFunction( "_mul",      new Func__mul<            Integer,           Real,       Real >() );
        addFunction( "_mul",      new Func__mul<               Real,        Integer,       Real >() );
        addFunction( "_mul",      new Func__mul<         MatrixReal,     MatrixReal, MatrixReal >() );
        addFunction( "_mul",      new Func__mul<         MatrixReal,           Real, MatrixReal >() );
        addFunction( "_mul",      new Func__mul<               Real,     MatrixReal, MatrixReal >() );
        addFunction( "_sub",      new Func__sub<            Integer,        Integer,    Integer >() );
        addFunction( "_sub",      new Func__sub<               Real,           Real,       Real >() );
        addFunction( "_sub",      new Func__sub<            Integer,           Real,       Real >() );
        addFunction( "_sub",      new Func__sub<               Real,        Integer,       Real >() );
        addFunction( "_sub",      new Func__sub<         MatrixReal,     MatrixReal, MatrixReal >() );
        
        /* Add basic logic templated functions */
        addFunction( "_and",      new Func__and<            Integer,        Integer >()             );
        addFunction( "_and",      new Func__and<               Real,           Real >()             );
        addFunction( "_and",      new Func__and<            Integer,           Real >()             );
        addFunction( "_and",      new Func__and<               Real,        Integer >()             );
        addFunction( "_and",      new Func__and<            Boolean,        Boolean >()             );
        addFunction( "_and",      new Func__and<            Boolean,        Integer >()             );
        addFunction( "_and",      new Func__and<            Boolean,           Real >()             );
        addFunction( "_and",      new Func__and<            Integer,        Boolean >()             );
        addFunction( "_and",      new Func__and<               Real,        Boolean >()             );
        addFunction( "_eq",       new Func__eq<             Integer,        Integer >()             );
        addFunction( "_eq",       new Func__eq<                Real,           Real >()             );
        addFunction( "_eq",       new Func__eq<             Integer,           Real >()             );
        addFunction( "_eq",       new Func__eq<                Real,        Integer >()             );
        addFunction( "_eq",       new Func__eq<             Boolean,        Boolean >()             );
        addFunction( "_eq",       new Func__eq<             Boolean,        Integer >()             );
        addFunction( "_eq",       new Func__eq<             Boolean,           Real >()             );
        addFunction( "_eq",       new Func__eq<             Integer,        Boolean >()             );
        addFunction( "_eq",       new Func__eq<                Real,        Boolean >()             );
        addFunction( "_ge",       new Func__ge<             Integer,        Integer >()             );
        addFunction( "_ge",       new Func__ge<                Real,           Real >()             );
        addFunction( "_ge",       new Func__ge<             Integer,           Real >()             );
        addFunction( "_ge",       new Func__ge<                Real,        Integer >()             );
        addFunction( "_ge",       new Func__ge<             Boolean,        Boolean >()             );
        addFunction( "_ge",       new Func__ge<             Boolean,        Integer >()             );
        addFunction( "_ge",       new Func__ge<             Boolean,           Real >()             );
        addFunction( "_ge",       new Func__ge<             Integer,        Boolean >()             );
        addFunction( "_ge",       new Func__ge<                Real,        Boolean >()             );
        addFunction( "_gt",       new Func__gt<             Integer,        Integer >()             );
        addFunction( "_gt",       new Func__gt<                Real,           Real >()             );
        addFunction( "_gt",       new Func__gt<             Integer,           Real >()             );
        addFunction( "_gt",       new Func__gt<                Real,        Integer >()             );
        addFunction( "_gt",       new Func__gt<             Boolean,        Boolean >()             );
        addFunction( "_gt",       new Func__gt<             Boolean,        Integer >()             );
        addFunction( "_gt",       new Func__gt<             Boolean,           Real >()             );
        addFunction( "_gt",       new Func__gt<             Integer,        Boolean >()             );
        addFunction( "_gt",       new Func__gt<                Real,        Boolean >()             );
        addFunction( "_lt",       new Func__lt<             Integer,        Integer >()             );
        addFunction( "_lt",       new Func__lt<                Real,           Real >()             );
        addFunction( "_lt",       new Func__lt<             Integer,           Real >()             );
        addFunction( "_lt",       new Func__lt<                Real,        Integer >()             );
        addFunction( "_lt",       new Func__lt<             Boolean,        Boolean >()             );
        addFunction( "_lt",       new Func__lt<             Boolean,        Integer >()             );
        addFunction( "_lt",       new Func__lt<             Boolean,           Real >()             );
        addFunction( "_lt",       new Func__lt<             Integer,        Boolean >()             );
        addFunction( "_lt",       new Func__lt<                Real,        Boolean >()             );
        addFunction( "_le",       new Func__le<             Integer,        Integer >()             );
        addFunction( "_le",       new Func__le<                Real,           Real >()             );
        addFunction( "_le",       new Func__le<             Integer,           Real >()             );
        addFunction( "_le",       new Func__le<                Real,        Integer >()             );
        addFunction( "_le",       new Func__le<             Boolean,        Boolean >()             );
        addFunction( "_le",       new Func__le<             Boolean,        Integer >()             );
        addFunction( "_le",       new Func__le<             Boolean,           Real >()             );
        addFunction( "_le",       new Func__le<             Integer,        Boolean >()             );
        addFunction( "_le",       new Func__le<                Real,        Boolean >()             );
        addFunction( "_ne",       new Func__ne<             Integer,        Integer >()             );
        addFunction( "_ne",       new Func__ne<                Real,           Real >()             );
        addFunction( "_ne",       new Func__ne<             Integer,           Real >()             );
        addFunction( "_ne",       new Func__ne<                Real,        Integer >()             );
        addFunction( "_ne",       new Func__ne<             Boolean,        Boolean >()             );
        addFunction( "_ne",       new Func__ne<             Boolean,        Integer >()             );
        addFunction( "_ne",       new Func__ne<             Boolean,           Real >()             );
        addFunction( "_ne",       new Func__ne<             Integer,        Boolean >()             );
        addFunction( "_ne",       new Func__ne<                Real,        Boolean >()             );
        addFunction( "_or",       new Func__or<             Integer,        Integer >()             );
        addFunction( "_or",       new Func__or<                Real,           Real >()             );
        addFunction( "_or",       new Func__or<             Integer,           Real >()             );
        addFunction( "_or",       new Func__or<                Real,        Integer >()             );
        addFunction( "_or",       new Func__or<             Boolean,        Boolean >()             );
        addFunction( "_or",       new Func__or<             Boolean,        Integer >()             );
        addFunction( "_or",       new Func__or<             Boolean,           Real >()             );
        addFunction( "_or",       new Func__or<             Integer,        Boolean >()             );
        addFunction( "_or",       new Func__or<                Real,        Boolean >()             );
        
        /* Add regular functions (alphabetical order) */
        addFunction( "clamp",     new Func_clamp()     ); 
        addFunction( "ls",        new Func_ls()        );
        addFunction( "model",     new Func_model()     );
        addFunction( "normalize", new Func_normalize() );
        addFunction( "q",         new Func_quit()      );
        addFunction( "quit",      new Func_quit()      );
        addFunction( "str",       new Func_str()       );
        
        /* Add math functions (alphabetical order) */ 
        addFunction( "expf",      new Func_exp()   );
        addFunction( "ln",        new Func_ln()    );
        addFunction( "log",       new Func_log()   );
        addFunction( "power",     new Func_power() );
        addFunction( "sin",       new Func_sin()   );
        addFunction( "sqrt",      new Func_sqrt()  );
        
        /* Add phylogeny-related functions */
        addFunction( "gtr",       new Func_gtr()          );
        addFunction( "read",      new Func_readAlignment() );

        /* Add regular templated functions (alphabetic order) */
        // TODO: Add v and s functions
        addFunction( "transpose", new Func_transpose< MatrixReal >() );
        addFunction( "v",         new Func_v< Real,     VectorReal    >() );
        addFunction( "v",         new Func_v< Integer,  VectorInteger >() );

    }
    catch(RbException& rbException) {

        RBOUT("Caught an exception while initializing the workspace\n");
        std::ostringstream msg;
        rbException.printValue(msg);
        msg << std::endl;
        RBOUT(msg.str());

        RBOUT("Please report this bug to the RevBayes Development Core Team");

        RBOUT("Press any character to exit the program.");
        getchar();
        exit(1);
    }
}


