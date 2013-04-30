#include "BranchHeterogeneousCharEvoModel.h"
#include "BetaDistribution.h"
#include "BetaSimplexMove.h"
#include "Clade.h"
#include "ConstantNode.h"
#include "ConstantBirthDeathProcess.h"
#include "ContinuousStochasticNode.h"
#include "DeterministicNode.h"
#include "DirichletDistribution.h"
#include "ExponentialDistribution.h"
#include "FileMonitor.h"
#include "FixedNodeheightPruneRegraft.h"
#include "GtrRateMatrixFunction.h"
#include "Mcmc.h"
#include "Model.h"
#include "Monitor.h"
#include "Move.h"
#include "NarrowExchange.h"
#include "NclReader.h"
#include "NearestNeighborInterchange.h"
#include "NodeTimeSlideBeta.h"
#include "NodeTimeSlideUniform.h"
#include "NucleotideFrequenciesFromGcContentFunction.h"
#include "RandomNumberGenerator.h"
#include "RateMatrix_GTR.h"
#include "RbFileManager.h"
#include "RootTimeSlide.h"
#include "ScaleMove.h"
#include "ScreenMonitor.h"
#include "SimplexMove.h"
#include "SlidingMove.h"
#include "SubtreeScale.h"
#include "TestBranchHeterogeneousTamura92Model.h"
#include "TimeTree.h"
#include "TreeHeightStatistic.h"
#include "TreeScale.h"
#include "UniformDistribution.h"
#include "Tamura92RateMatrixFunction.h"
#include "VectorFunction.h"
#include "RbVectorFunction.h"

using namespace RevBayesCore;

TestBranchHeterogeneousTamura92Model::TestBranchHeterogeneousTamura92Model(const std::string &afn, const std::string &tFn, int gen) : alignmentFilename( afn ), treeFilename( tFn ), mcmcGenerations( gen ){
    
}

TestBranchHeterogeneousTamura92Model::~TestBranchHeterogeneousTamura92Model() {
    // nothing to do
}



bool TestBranchHeterogeneousTamura92Model::run( void ) {
    
    std::vector<unsigned int> seeds;
    seeds.push_back(7);
    seeds.push_back(4);
    GLOBAL_RNG->setSeed( seeds );
    
    /* First, we read in the data */
    // the matrix
    std::vector<AbstractCharacterData*> data = NclReader::getInstance().readMatrices(alignmentFilename);
    std::cout << "Read " << data.size() << " matrices." << std::endl;
    std::cout << data[0] << std::endl;
    
    std::vector<TimeTree*> trees = NclReader::getInstance().readTimeTrees( treeFilename );
    std::cout << "Read " << trees.size() << " trees." << std::endl;
    std::cout << trees[0]->getNewickRepresentation() << std::endl;
    
    
    /* set up the model graph */
    
    //////////////////////
    // first the priors //
    //////////////////////
    
    // birth-death process priors
    StochasticNode<double> *div = new StochasticNode<double>("diversification", new UniformDistribution(new ConstantNode<double>("div_lower", new double(0.0)), new ConstantNode<double>("div_upper", new double(100.0)) ));
    ConstantNode<double> *turn = new ConstantNode<double>("turnover", new double(0.0));
    ConstantNode<double> *rho = new ConstantNode<double>("rho", new double(1.0));
    
    //ts/tv ratio:
    ConstantNode<double > *tstv_prior = new ConstantNode<double >( "tstv_prior", new double(0.25) );
    ContinuousStochasticNode *tstv = new ContinuousStochasticNode("tstv", new ExponentialDistribution(tstv_prior) );

    //GC content prior:
    ConstantNode<double > *eq_gc_prior = new ConstantNode<double >( "eq_gc_prior_ab", new double(1.0) );    

    //Equilibrium GC frequency: one per branch, see below.
    //ContinuousStochasticNode *eqgc = new ContinuousStochasticNode("eq_gc", new BetaDistribution(eq_gc_prior,eq_gc_prior) );
    
    //Root GC frequency
    
    StochasticNode< double  > *omega = new StochasticNode< double >( "omega", new BetaDistribution(eq_gc_prior,eq_gc_prior) );
    DeterministicNode<std::vector<double> > *rf = new DeterministicNode< std::vector<double> >( "rf", new NucleotideFrequenciesFromGcContentFunction( omega ) );

    std::cout << "omega:\t" << omega->getValue() << std::endl;
    std::cout << "rf:\t" << rf->getValue() << std::endl;
    std::cout << "tstv:\t" << tstv->getValue() << std::endl;
    
    //Declaring a vector of GTR matrices
    size_t numBranches = 2*data[0]->getNumberOfTaxa() - 2;
    std::vector<ContinuousStochasticNode*> thetas;
    //  std::vector<StochasticNode < std::vector<double> >* > ers;
    std::vector< const TypedDagNode < RateMatrix >* > qs;
    
    for (unsigned int i = 0 ; i < numBranches ; i++ ) {
        std::ostringstream eq_gc_name;
        eq_gc_name << "eq_gc(" << i << ")";
        thetas.push_back(new ContinuousStochasticNode( eq_gc_name.str(), new BetaDistribution(eq_gc_prior,eq_gc_prior) ) );
        
        //std::cout << "theta " << i << ":\t" << thetas[i]->getValue() << "\n"; 

        //  ers.push_back(new StochasticNode<std::vector<double> >( "er", new DirichletDistribution(e) ) );
        std::ostringstream q_name;
        q_name << "q(" << i << ")";

        
        qs.push_back(new DeterministicNode<RateMatrix>( q_name.str(), new Tamura92RateMatrixFunction( thetas[i], tstv) ));

        //std::cout << "Q:\t" << qs[i]->getValue() << std::endl;
    }
    
    DeterministicNode< RbVector< RateMatrix > >* qs_node = new DeterministicNode< RbVector< RateMatrix > >( "q_vector", new RbVectorFunction<RateMatrix>(qs) );
    
    std::vector<std::string> names = data[0]->getTaxonNames();
    ConstantNode<std::vector<double> > *met = new ConstantNode<std::vector<double> >("MET",new std::vector<double>() );
    ConstantNode<std::vector<double> > *mep = new ConstantNode<std::vector<double> >("MESP",new std::vector<double>() );
    StochasticNode<TimeTree> *tau = new StochasticNode<TimeTree>( "tau", new ConstantBirthDeathProcess(div, turn, met, mep, rho, "uniform", "survival", int(names.size()), names, std::vector<Clade>()) );
    
    tau->setValue( trees[0] );
    std::cout << "tau:\t" << tau->getValue() << std::endl;
    
    // and the character model
   StochasticNode< AbstractCharacterData > *charactermodel = new StochasticNode< AbstractCharacterData >("S", new BranchHeterogeneousCharEvoModel<DnaState, TimeTree>(tau, qs_node, rf, true, data[0]->getNumberOfCharacters()) );
   charactermodel->clamp( data[0] );
    
    std::cout <<"Hehe "<<std::endl;

    
    /* add the moves */
    std::vector<Move*> moves;
    moves.push_back( new ScaleMove(div, 1.0, true, 2.0) );
    moves.push_back( new NearestNeighborInterchange( tau, 5.0 ) );
    moves.push_back( new NarrowExchange( tau, 10.0 ) );
    moves.push_back( new FixedNodeheightPruneRegraft( tau, 2.0 ) );
    moves.push_back( new SubtreeScale( tau, 5.0 ) );
    moves.push_back( new TreeScale( tau, 1.0, true, 2.0 ) );
    moves.push_back( new NodeTimeSlideUniform( tau, 30.0 ) );
    moves.push_back( new RootTimeSlide( tau, 1.0, true, 2.0 ) );
//    moves.push_back( new SimplexMove( er, 10.0, 1, true, 2.0 ) );
//    moves.push_back( new SimplexMove( er, 100.0, 6, true, 2.0 ) );
    moves.push_back( new BetaSimplexMove( omega, 10.0, true, 2.0 ) );
    moves.push_back( new ScaleMove(tstv, 1.0, true, 2.0) );

    for (unsigned int i = 0 ; i < numBranches ; i ++ ) {
        moves.push_back( new BetaSimplexMove( dynamic_cast<StochasticNode<double>* >(thetas[i]), 10.0, true, 2.0 ) );
        moves.push_back( new SlidingMove( thetas[i], 0.05, true, 2.0) );
        //    moves.push_back( new SimplexMove( ers[i], 100.0, 6, true, 2.0 ) );
//        moves.push_back( new SimplexMove( pis[i], 100.0, 4, true, 2.0 ) );
    }
    
    // add some tree stats to monitor
    DeterministicNode<double> *treeHeight = new DeterministicNode<double>("TreeHeight", new TreeHeightStatistic(tau) );
    
    /* add the monitors */
    std::vector<Monitor*> monitors;
    std::set<DagNode*> monitoredNodes;
    //    monitoredNodes.insert( er );
    //    monitoredNodes.insert( pi );
    monitoredNodes.insert( div );
    monitors.push_back( new FileMonitor( monitoredNodes, 10, "TestBranchHeterogeneousTamura92Model.log", "\t" ) );
    std::set<DagNode*> monitoredNodes1;
//    monitoredNodes1.insert( er );
    for (unsigned int i = 0 ; i < numBranches ; i ++ ) {
        monitoredNodes1.insert( thetas[i] );
    }
    monitoredNodes1.insert( rf );
    monitoredNodes1.insert( treeHeight );
    monitors.push_back( new FileMonitor( monitoredNodes1, 10, "TestBranchHeterogeneousTamura92ModelSubstRates.log", "\t" ) );
    monitors.push_back( new ScreenMonitor( monitoredNodes1, 10, "\t" ) );
    std::set<DagNode*> monitoredNodes2;
    monitoredNodes2.insert( tau );
    monitors.push_back( new FileMonitor( monitoredNodes2, 10, "TestBranchHeterogeneousTamura92Model.tree", "\t", false, false, false ) );
    
    /* instantiate the model */
    Model myModel = Model(qs[0]);
    
    /* instiate and run the MCMC */
    Mcmc myMcmc = Mcmc( myModel, moves, monitors );
    myMcmc.run(mcmcGenerations);
    
    myMcmc.printOperatorSummary();
    
    /* clean up */
    //    for (size_t i = 0; i < 10; ++i) {
    //        delete x[i];
    //    }
    //    delete [] x;
    delete div;
    //    delete sigma;
    //    delete a;
    //    delete b;
    //    delete c;
    for (std::vector<Move*>::iterator it = moves.begin(); it != moves.end(); ++it) {
        const Move *theMove = *it;
        delete theMove;
    }
    for (std::vector<Monitor*>::iterator it = monitors.begin(); it != monitors.end(); ++it) {
        const Monitor *theMonitor = *it;
        delete theMonitor;
    }
    
    std::cout << "Finished GTR model test." << std::endl;
    
    return true;
}
