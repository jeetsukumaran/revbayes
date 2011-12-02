/**
 * @file
 * This file contains the declaration of NclReader, which is the wrapper class for the NCL 
 * for reading character matrices, sequences, and trees.
 
 * @brief Declaration of NclReader
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2009-11-19 17:29:33 +0100 (Tor, 19 Nov 2009) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-08-16, version 1.0
 * @extends DAGNode
 *
 * $Id$
 */

#ifndef NclReader_H
#define NclReader_H

#include "ncl.h"
#include "nxsmultiformat.h"
#include "RbPtr.h"
#include "TopologyNode.h"
#include "TreePlate.h"
#include <map>
#include <set>
#include <string>
#include <vector>

class CharacterData;
class Topology;

class NclReader{
    
        friend class NxsBlock;
        friend class NxsUnalignedBlock;
        
    public:
        void                                    addWarning(std::string s) { warningsSummary.insert(s); }                        //!< Add a warning to the warnings vector
        void                                    clearWarnings(void) { warningsSummary.clear(); }                                //!< Clear all of the warnings from the warnings vector
        size_t                                  getNumWarnings(void) { return warningsSummary.size(); }                         //!< Return the number of warnings
        std::set<std::string>&                  getWarnings(void) { return warningsSummary; }                                   //!< Get a reference to the warnings vector
        static NclReader&                       getInstance(void);                                                              //!< Get a reference to this singleton class
        std::vector<RbPtr<CharacterData> >      readMatrices(const std::map<std::string,std::string>& fileMap);                 //!< Read a list of file names contained in a map (with file format info too)
        std::vector<RbPtr<CharacterData> >      readMatrices(const std::vector<std::string> fn, const std::string fileFormat, const std::string dataType, const bool isInterleaved);              //!< Read a list of file names contained in a vector of strings
        
        bool                                    isFastaFile(std::string& fn, std::string& dType);                               //!< Checks if the file is in Fasta format
        bool                                    isNexusFile(std::string& fn);                                                   //!< Checks if the file is in NEXUS format
        bool                                    isPhylipFile(std::string& fn, std::string& dType, bool& isInterleaved);         //!< Checks if the file is in Phylip format

        // TAH: stuff for reading trees
        RbPtr<std::vector<RbPtr<TreePlate> > >  readTrees(const std::string fn, const std::string fileFormat);                  //!< Read trees
        void                                    clearContent(void) { nexusReader.ClearContent(); }                              //!< Clear the content of the NCL object
        
    private:
                                                NclReader(void) { }                                                             //!< Default constructor
                                                NclReader(const NclReader& r) { }                                               //!< Copy constructor
        virtual                                ~NclReader(void) { }                                                             //!< Destructor
        
        RbPtr<CharacterData>                    createAminoAcidMatrix(NxsCharactersBlock* charblock);                           //!< Create an object to hold amino acid data
        RbPtr<CharacterData>                    createContinuousMatrix(NxsCharactersBlock* charblock);                          //!< Create an object to hold continuous data
        RbPtr<CharacterData>                    createDnaMatrix(NxsCharactersBlock* charblock);                                 //!< Create an object to hold DNA data
        RbPtr<CharacterData>                    createRnaMatrix(NxsCharactersBlock* charblock);                                 //!< Create an object to hold RNA data
        RbPtr<CharacterData>                    createStandardMatrix(NxsCharactersBlock* charblock);                            //!< Create an object to hold standard data
        RbPtr<CharacterData>                    createUnalignedAminoAcidMatrix(NxsUnalignedBlock* charblock);                   //!< Create an object to hold amino acid data
        RbPtr<CharacterData>                    createUnalignedDnaMatrix(NxsUnalignedBlock* charblock);                         //!< Create an object to hold DNA data
        bool                                    fileExists(const char *fn) const;                                               //!< Returns whether a file exists
        std::string                             findFileNameFromPath(const std::string& fp) const;                              //!< Returns the file name from a file path
        std::string                             intuitDataType(std::string& s);                                                 //!< Attempt to determine the type of data

        // methods for reading sequence alignments
        std::vector<RbPtr<CharacterData> >      convertFromNcl(std::vector<std::string>& fnv);                                  //!< Reads the blocks stored by NCL and converts them to RevBayes character matrices 
        std::vector<RbPtr<CharacterData> >      readMatrices(const char* fileName, const std::string fileFormat, const std::string dataType, const bool isInterleaved);                          //!< Reads a single file using NCL
        void                                    setExcluded(const NxsCharactersBlock* charblock, RbPtr<CharacterData> cMat ) const;       //!< Set excluded taxa and excluded characters

        // methods for reading trees
        void                                    constructTreefromNclRecursively(RbPtr<TopologyNode> tn, const NxsSimpleNode* tnNcl);  //!< Constructs a tree from NCL
        RbPtr<std::vector<RbPtr<TreePlate> > >  readTrees(const char* fileName, const std::string fileFormat);                  //!< Reads trees contained in a file
        RbPtr<std::vector<RbPtr<TreePlate> > >  convertTreesFromNcl(void);                                                      //!< Converts trees stored by NCL into RevBayes formatted trees
        RbPtr<TreePlate>                        translateNclSimpleTreeToTree(NxsSimpleTree &nTree);                             //!< Translate a single NCL tree into a RevBayes tree
        
        MultiFormatReader                       nexusReader;                                                                    //!< The NCL object that reads the files
        std::set<std::string>                   warningsSummary;                                                                //!< A vector that contains the warnings that acumulate
};

#endif