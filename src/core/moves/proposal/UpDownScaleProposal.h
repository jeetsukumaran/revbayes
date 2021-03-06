#ifndef UpDownScaleProposal_H
#define UpDownScaleProposal_H

#include "RbVector.h"
#include "Proposal.h"
#include "StochasticNode.h"
#include "TimeTree.h"

#include <ostream>
#include <vector>
#include <string>

namespace RevBayesCore {
    
    /**
     * @brief Up-Down Scaling Proposal for several parameters jointly.
     *
     *
     * This proposal randomly scales all a set of parameter up while the other set of parameters
     * is scaled down by the same value. This should hopefully improve mixing in many cases.
     * The actual scaling factor is computed by sf = exp( lambda * ( u - 0.5 ) )
     * where u ~ Uniform(0,1).
     *
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @copyright GPL version 3
     * @since 2015-05-25, version 1.0
     *
     */
    class UpDownScaleProposal : public Proposal {
        
    public:
        UpDownScaleProposal(double l);                                 //!< Constructor
        
        void                                        addDownVariable(StochasticNode<double> *v);                                         //!< Add an up-scaling variable
        void                                        addDownVariable(StochasticNode<RbVector<double> > *v);                              //!< Add an up-scaling variable
        void                                        addDownVariable(StochasticNode<TimeTree> *v);                                       //!< Add an up-scaling variable
        void                                        addUpVariable(StochasticNode<double> *v);                                           //!< Add an up-scaling variable
        void                                        addUpVariable(StochasticNode<RbVector<double> > *v);                                //!< Add an up-scaling variable
        void                                        addUpVariable(StochasticNode<TimeTree> *v);                                         //!< Add an up-scaling variable
        void                                        cleanProposal(void);                                                                //!< Clean up proposal
        UpDownScaleProposal*                        clone(void) const;                                                                  //!< Clone object
        double                                      doProposal(void);                                                                   //!< Perform proposal
        const std::string&                          getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        void                                        printParameterSummary(std::ostream &o) const;                                       //!< Print the parameter summary
        void                                        prepareProposal(void);                                                              //!< Prepare the proposal
        void                                        tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                        undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                        swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        // parameters
        
        std::vector<StochasticNode<double> *>               upScalarVariables;
        std::vector<StochasticNode<RbVector<double> > *>    upVectorVariables;
        std::vector<StochasticNode<TimeTree> *>             upTreeVariables;
        std::vector<StochasticNode<double> *>               downScalarVariables;
        std::vector<StochasticNode<RbVector<double> > *>    downVectorVariables;
        std::vector<StochasticNode<TimeTree> *>             downTreeVariables;
        
        double                                              lambda;                                                                             //!< The scale parameter of the Proposal (larger lambda -> larger proposals).
        double                                              storedScalingFactor;                                                                        //!< The stored value of the last modified UpDown.
        
        
        
    };
    
}

#endif

