#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "Ellipsis.h"
#include "Func_readAncestralStateTrace.h"
#include "ModelVector.h"
#include "NclReader.h"
#include "NewickConverter.h"
#include "OptionRule.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RlString.h"
#include "RlAncestralStateTrace.h"
#include "RlUtils.h"
#include "StringUtilities.h"
#include "AncestralStateTrace.h"
#include "TreeUtilities.h"
#include "UserInterface.h"

#include <map>
#include <set>
#include <sstream>


using namespace RevLanguage;

/** Clone object */
Func_readAncestralStateTrace* Func_readAncestralStateTrace::clone( void ) const
{
    
    return new Func_readAncestralStateTrace( *this );
}


/** Execute function */
RevPtr<Variable> Func_readAncestralStateTrace::execute( void ) {
    
    // get the information from the arguments for reading the file
    const std::string&  fn       = static_cast<const RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    const std::string&  sep      = static_cast<const RlString&>( args[1].getVariable()->getRevObject() ).getValue();
    
    // check that the file/path name has been correctly specified
    RevBayesCore::RbFileManager myFileManager( fn );
    
    if ( !myFileManager.testFile() || !myFileManager.testDirectory() )
    {
        std::string errorStr = "";
        myFileManager.formatError(errorStr);
        throw RbException(errorStr);
    }
    
    // set up a vector of strings containing the name or names of the files to be read
    std::vector<std::string> vectorOfFileNames;
    if ( myFileManager.isFile() )
    {
        vectorOfFileNames.push_back( myFileManager.getFullFileName() );
    }
    else
    {
        myFileManager.setStringWithNamesOfFilesInDirectory( vectorOfFileNames );
    }
    
    RevObject *rv;
	rv = readAncestralStates(vectorOfFileNames, sep);
    
    
    return new Variable( rv );
}


/** Format the error exception string for problems specifying the file/path name */
void Func_readAncestralStateTrace::formatError(RevBayesCore::RbFileManager& fm, std::string& errorStr)
{
    
    bool fileNameProvided    = fm.isFileNamePresent();
    bool isFileNameGood      = fm.testFile();
    bool isDirectoryNameGood = fm.testDirectory();
    
    if ( fileNameProvided == false && isDirectoryNameGood == false )
    {
        errorStr += "Could not read contents of directory \"" + fm.getFilePath() + "\" because the directory does not exist";
    }
    else if (fileNameProvided == true && (isFileNameGood == false || isDirectoryNameGood == false)) {
        errorStr += "Could not read file named \"" + fm.getFileName() + "\" in directory named \"" + fm.getFilePath() + "\" ";
        if (isFileNameGood == false && isDirectoryNameGood == true)
            errorStr += "because the file does not exist";
        else if (isFileNameGood == true && isDirectoryNameGood == false)
            errorStr += "because the directory does not exist";
        else
            errorStr += "because neither the directory nor the file exist";
    }
}


/** Get argument rules */
const ArgumentRules& Func_readAncestralStateTrace::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rulesSet = false;
    
    if (!rulesSet)
    {
		
        argumentRules.push_back( new ArgumentRule( "file"     , RlString::getClassTypeSpec(), ArgumentRule::BY_VALUE ) );
        argumentRules.push_back( new ArgumentRule( "separator", RlString::getClassTypeSpec(), ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("\t") ) );
        rulesSet = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readAncestralStateTrace::getClassType(void) {
    
    static std::string revType = "Func_readAncestralStateTrace";
    
	return revType;
}

/** Get class type spec describing type of object */
const TypeSpec& Func_readAncestralStateTrace::getClassTypeSpec(void) {
    
    static TypeSpec revTypeSpec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return revTypeSpec;
}

/** Get type spec */
const TypeSpec& Func_readAncestralStateTrace::getTypeSpec( void ) const {
    
    static TypeSpec typeSpec = getClassTypeSpec();
    
    return typeSpec;
}


/** Get return type */
const TypeSpec& Func_readAncestralStateTrace::getReturnType( void ) const {
    
    static TypeSpec returnTypeSpec = RevLanguage::AncestralStateTrace::getClassTypeSpec();
    return returnTypeSpec;
}

// convert vector of strings to RL Variable
//        const std::vector<std::string>& n = this->dagNode->getValue().getTipNames();
//		  return new Variable( new ModelVector<RlString>( n ) );


AncestralStateTrace* Func_readAncestralStateTrace::readAncestralStates(const std::vector<std::string> &vectorOfFileNames, const std::string &delimitter)
{
    
    
    std::vector<RevBayesCore::AncestralStateTrace> data;
    
    
    // Set up a map with the file name to be read as the key and the file type as the value. Note that we may not
    // read all of the files in the string called "vectorOfFileNames" because some of them may not be in a format
    // that can be read.
    std::map<std::string,std::string> fileMap;
    for (std::vector<std::string>::const_iterator p = vectorOfFileNames.begin(); p != vectorOfFileNames.end(); p++)
    {
        bool hasHeaderBeenRead = false;
        const std::string &fn = *p;
        
        /* Open file */
        std::ifstream inFile( fn.c_str() );
        
        if ( !inFile )
            throw RbException( "Could not open file \"" + fn + "\"" );
        
        /* Initialize */
        std::string commandLine;
        std::cerr << "Processing file \"" << fn << "\"" << std::endl;
        
        size_t index = 0;
        
        /* Command-processing loop */
        while ( inFile.good() )
        {
            
            // Read a line
            std::string line;
            getline( inFile, line );
            
            // skip empty lines
            //line = StringByTrimmingCharactersInSet:[NSCharacterSet whitespaceAndNewlineCharacterSet]];
            if (line.length() == 0)
            {
                continue;
            }
            
            
            // removing comments
            if (line[0] == '#') {
                continue;
            }
            
            // splitting every line into its columns
            std::vector<std::string> columns;
            
            // we should provide other delimiters too
            StringUtilities::stringSplit(line, delimitter, columns);
            
            // we assume a header at the first line of the file
            if (!hasHeaderBeenRead) {
                
                for (size_t j=1; j<columns.size(); j++) {
                    
                    std::string parmName = columns[j];
                    if ( parmName == "Posterior" || parmName == "Likelihood" || parmName == "Prior") {
                        continue;
                    }
                    index = j;
                    
                    RevBayesCore::AncestralStateTrace t = RevBayesCore::AncestralStateTrace();
                    
                    t.setParameterName(parmName);
                    t.setFileName(fn);
                    
                    data.push_back( t );
                }
                
                hasHeaderBeenRead = true;
                
                continue;
            }
            
            // adding values to the Traces
            RevBayesCore::AncestralStateTrace& t = data[0];
            
            std::string anc_state = columns[index];
            
            t.addObject( anc_state );
        }
    }
    
    return new RevLanguage::AncestralStateTrace( data[0] );
}




