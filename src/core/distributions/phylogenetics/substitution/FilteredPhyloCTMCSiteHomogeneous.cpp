#include <cstdlib>
#include <vector>
#include <cmath>
#include <FilteredPhyloCTMCSiteHomogeneous.h>
using namespace RevBayesCore;
class RevBayesCore::AscertainmentBiasCorrectionStruct {
    public:
        void computeInternalAscBias(const AscertainmentBiasCorrectionStruct * ascLeft, const AscertainmentBiasCorrectionStruct * ascRight, size_t numSiteRates, size_t numStates, size_t numPatterns, const double ** tpMats);
        void computeTipAscBias(size_t numSiteRates, size_t numStates, size_t numPatterns, const double ** tpMats,  const std::vector<bool> &gap_node, const std::vector<unsigned long> &char_node, bool usingAmbiguousCharacters);
        double computeAscBiasLnProbCorrection2Node(const AscertainmentBiasCorrectionStruct * ascRight, const size_t numSiteRates, const double *rootFreq, const size_t numStates, const size_t * patternCounts, const double p_inv,const std::vector<bool> &  siteInvariant, const std::vector<size_t> & invariantSiteIndex) const;
        double computeAscBiasLnProbCorrection3Node(const AscertainmentBiasCorrectionStruct * ascRight, const AscertainmentBiasCorrectionStruct * ascMiddle, const size_t numSiteRates, const double *rootFreq, const size_t numStates, const size_t * patternCounts, const double p_inv,const std::vector<bool> &  siteInvariant, const std::vector<size_t> & invariantSiteIndex) const;
};

std::vector<RevBayesCore::AscertainmentBiasCorrectionStruct *> RevBayesCore::allocAscBiasCorrStructs(const size_t numCopies, const size_t numNodes) {
    std::vector<RevBayesCore::AscertainmentBiasCorrectionStruct *> x(numCopies*numNodes, 0L);
    try{ 
        for (size_t i = 0 ; i < numCopies*numNodes; ++i) {
            x[i] = new RevBayesCore::AscertainmentBiasCorrectionStruct();
        }
    } catch (...) {
        freeAscBiasCorrStructs(x);
        throw;
    }
}

void RevBayesCore::freeAscBiasCorrStructs(std::vector<RevBayesCore::AscertainmentBiasCorrectionStruct *> &x) {
    for (size_t i = 0 ; i < x.size(); ++i) {
        delete x[i];
        x[i] = NULL;
    }
    x.clear();
}
double RevBayesCore::computeRootFilteredLikelihood2Nodes(const double *p_left,
                                   const double *p_right,
                                   const size_t numSiteRates,
                                   const double * rootFreq,
                                   const size_t numStates,
                                   const size_t * patternCounts,
                                   const size_t numPatterns,
                                   const size_t siteOffset,
                                   const size_t mixtureOffset,
                                   const double p_inv,
                                   const std::vector<bool> & siteInvariant,
                                   const std::vector<size_t> & invariantSiteIndex,
                                   const AscertainmentBiasCorrectionStruct *ascLeft,
                                   const AscertainmentBiasCorrectionStruct *ascRight,
                                   double & uncorrectedLnProb,
                                   double & ascBiasLnProb
                                   ) {
    uncorrectedLnProb = computeRootLikelihood2Nodes(p_left,
                                                      p_right,
                                                      numSiteRates,
                                                      rootFreq,
                                                      numStates,
                                                      patternCounts,
                                                      numPatterns,
                                                      siteOffset,
                                                      mixtureOffset,
                                                      p_inv,
                                                      siteInvariant,
                                                      invariantSiteIndex);
    ascBiasLnProb = ascLeft->computeAscBiasLnProbCorrection2Node(ascRight, numSiteRates, rootFreq, numStates, patternCounts, p_inv, siteInvariant, invariantSiteIndex);
    return ascBiasLnProb + uncorrectedLnProb;
}
double RevBayesCore::computeRootFilteredLikelihood3Nodes(const double *p_left,
                                   const double *p_right,
                                   const double *p_middle,
                                   const size_t numSiteRates,
                                   const double * rootFreq,
                                   const size_t numStates,
                                   const size_t * patternCounts,
                                   const size_t numPatterns,
                                   const size_t siteOffset,
                                   const size_t mixtureOffset,
                                   const double p_inv,
                                   const std::vector<bool> & siteInvariant,
                                   const std::vector<size_t> & invariantSiteIndex,
                                   const AscertainmentBiasCorrectionStruct *ascLeft,
                                   const AscertainmentBiasCorrectionStruct *ascRight,
                                   const AscertainmentBiasCorrectionStruct *ascMiddle,
                                   double & uncorrectedLnProb,
                                   double & ascBiasLnProb) {
    uncorrectedLnProb = RevBayesCore::computeRootLikelihood3Nodes(p_left,
                                                 p_right,
                                                 p_middle,
                                                 numSiteRates,
                                                 rootFreq,
                                                 numStates,
                                                 patternCounts,
                                                 numPatterns,
                                                 siteOffset,
                                                 mixtureOffset,
                                                 p_inv,
                                                 siteInvariant,
                                                 invariantSiteIndex);
    ascBiasLnProb = ascLeft->computeAscBiasLnProbCorrection3Node(ascRight, ascMiddle, numSiteRates, rootFreq, numStates, patternCounts, p_inv, siteInvariant, invariantSiteIndex);
    return ascBiasLnProb + uncorrectedLnProb;
}
void RevBayesCore::computeInternalNodeFilteredLikelihood(double * p_node,
                                            AscertainmentBiasCorrectionStruct *ascNode,
                                    const double *p_left,
                                    const double *p_right,
                                    const size_t numSiteRates,
                                    const size_t numStates,
                                    const size_t numPatterns,
                                    const size_t siteOffset,
                                    const size_t mixtureOffset,
                                    const double ** tpMats,
                                   const AscertainmentBiasCorrectionStruct *ascLeft,
                                   const AscertainmentBiasCorrectionStruct *ascRight) {
    computeInternalNodeLikelihood(p_node, p_left, p_right, numSiteRates, numStates, numPatterns, siteOffset, mixtureOffset, tpMats);
    ascNode->computeInternalAscBias(ascLeft, ascRight, numSiteRates, numStates, numPatterns, tpMats);
}
void RevBayesCore::computeTipNodeFilteredLikelihood(double * p_node,
                                            AscertainmentBiasCorrectionStruct *ascNode,
                               const size_t numSiteRates,
                               const size_t numStates,
                               const size_t numPatterns,
                               const size_t siteOffset,
                               const size_t mixtureOffset,
                               const double ** tpMats,
                               const std::vector<bool> &gap_node,
                               const std::vector<unsigned long> &char_node,
                               const bool usingAmbiguousCharacters) {
    computeTipNodeLikelihood(p_node, numSiteRates, numStates, numPatterns, siteOffset, mixtureOffset, tpMats, gap_node, char_node, usingAmbiguousCharacters);
    ascNode->computeTipAscBias(numSiteRates, numStates, numPatterns, tpMats, gap_node, char_node, usingAmbiguousCharacters);
}
