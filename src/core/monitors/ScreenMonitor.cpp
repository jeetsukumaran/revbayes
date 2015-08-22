/**
 * @file
 * This file contains the implementation of ScreenMonitor, used to save
 * information to the screen about the monitoring of a variable DAG node.
 *
 * @brief Implementation of ScreenMonitor
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-11-15 16:03:33 +0100 (Thu, 15 Nov 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-06-21, version 1.0
 *
 * $Id: ScreenMonitor.cpp 1833 2012-11-15 15:03:33Z hoehna $
 */


#include "ScreenMonitor.h"
#include "DagNode.h"
#include "Mcmc.h"
#include "Model.h"
#include "Monitor.h"
#include "RlUserInterface.h"
#include "StringUtilities.h"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <time.h>

using namespace RevBayesCore;

/* Constructor */
ScreenMonitor::ScreenMonitor(DagNode *n, int g, bool pp, bool l, bool pr) : Monitor(g,n),
    posterior( pp ),
    prior( pr ),
    likelihood( l ),
    printWaitingTime( true ),
    printElapsedTime( true ),
    printChainIndex( true ),
    prefixSeparator("   "),
    suffixSeparator("   |"),
    headerPrintingInterval( 20 ),
    startTime( NULL ),
    numCycles( 0 ),
    currentGen( 0 ),
    startGen( 0 ),
    replicateIndex( 0 ),
    enabled( true )
{
    
}


ScreenMonitor::ScreenMonitor(const std::vector<DagNode *> &n, int g, bool pp, bool l, bool pr, bool pci) : Monitor(g,n),
    posterior( pp ),
    prior( pr ),
    likelihood( l ),
    printWaitingTime( true ),
    printElapsedTime( true ),
//    printChainIndex( pci ),
    printChainIndex( true ),   // TODO modify revfunction to set this value
    prefixSeparator("   "),
    suffixSeparator("   |"),
    headerPrintingInterval( 20 ),
    startTime( NULL ),
    numCycles( 0 ),
    currentGen( 0 ),
    startGen( 0 ),
    replicateIndex( 0 ),
    enabled( true )
{
    
}


/* Clone the object */
ScreenMonitor* ScreenMonitor::clone(void) const
{
    
    return new ScreenMonitor(*this);
}



/** Monitor value at generation gen */
void ScreenMonitor::monitor(unsigned long gen)
{

    // only the enabled replacte and the cold chain is allowed to print
    if ( enabled == true  && mcmc->getChainPosteriorHeat() == 1.0 )
    {
        
        // get the printing frequency
        unsigned long samplingFrequency = printgen;
        
        if ( gen > 0 && gen % (headerPrintingInterval*samplingFrequency) == 0 )
        {
            printHeader();
        }
        
        if (gen % samplingFrequency == 0)
        {
            // handy string and stringstream
            std::string s;
            std::stringstream ss;
            
            // set column width
            int columnWidth = 12;
            
            // set cycle column width
            int cycleWidth = floor( log10( numCycles ) ) + 1;
            cycleWidth = 9 > cycleWidth ? 9 : cycleWidth;
            
            // final string to print for this iteration
            std::string output = "";

            // print the cycle number
            ss << gen;
            s = ss.str();
            StringUtilities::fillWithSpaces(s, cycleWidth, true);
            output = s + suffixSeparator;
            ss.str("");

            if ( printChainIndex )
            {
                ss << mcmc->getChainIndex();
                s = ss.str();
                StringUtilities::fillWithSpaces( s, columnWidth, false );
                output += prefixSeparator + s + suffixSeparator;
                ss.str("");
            }
            
            if ( posterior )
            {
                const std::vector<DagNode*> &n = model->getDagNodes();
                double pp = 0.0;
                for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it)
                {
                    pp += (*it)->getLnProbability();
                }
                
                ss << pp;
                s = ss.str();
                StringUtilities::fillWithSpaces( s, columnWidth, false );
                output += prefixSeparator + s + suffixSeparator;
                ss.str("");
            }
            
            if ( likelihood )
            {
                const std::vector<DagNode*> &n = model->getDagNodes();
                double pp = 0.0;
                for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it)
                {
                    if ( (*it)->isClamped() )
                    {
                        pp += (*it)->getLnProbability();
                    }
                }

                ss << pp;
                s = ss.str();
                StringUtilities::fillWithSpaces( s, columnWidth, false );
                output += prefixSeparator + s + suffixSeparator;
                ss.str("");
            }
            
            if ( prior )
            {
                const std::vector<DagNode*> &n = model->getDagNodes();
                double pp = 0.0;
                for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it)
                {
                    if ( !(*it)->isClamped() )
                    {
                        pp += (*it)->getLnProbability();
                    }
                }

                ss << pp;
                s = ss.str();
                StringUtilities::fillWithSpaces( s, columnWidth, false );
                output += prefixSeparator + s + suffixSeparator;
                ss.str("");
            }
            
            for (std::vector<DagNode*>::iterator i = nodes.begin(); i != nodes.end(); ++i)
            {
                // get the node
                DagNode *node = *i;
                
                // print the value
                node->printValueElements(ss, prefixSeparator + suffixSeparator, int( columnWidth ), false);
                output += prefixSeparator + ss.str() + suffixSeparator;
                ss.str("");
            }
            
            if ( printElapsedTime )
            {
                
                size_t timeUsed = time(NULL) - startTime;
                
                size_t hours   = timeUsed / 3600;
                size_t minutes = timeUsed / 60 - hours * 60;
                size_t seconds = timeUsed - minutes * 60 - hours * 3600;
                
                ss << std::setw( 2 ) << std::setfill( '0' ) << hours << ":";
                ss << std::setw( 2 ) << std::setfill( '0' ) << minutes << ":";
                ss << std::setw( 2 ) << std::setfill( '0' ) << seconds;
                
                output += prefixSeparator + ss.str() + suffixSeparator;
                ss.str("");
                
            }
            
            if ( printWaitingTime )
            {

                if ( (gen-startGen) <= samplingFrequency )
                {
                    output += prefixSeparator + "--:--:--" + suffixSeparator;
                }
                else
                {
                    double done = gen - startGen;
                    double total = numCycles - startGen;
                    double progress = done / total;
                    size_t timeUsed = time(NULL) - startTime;
                
                    size_t waitTime = double( timeUsed ) / progress - double( timeUsed );
                    
                    size_t hours   = waitTime / 3600;
                    size_t minutes = waitTime / 60 - hours * 60;
                    size_t seconds = waitTime - minutes * 60 - hours * 3600;
                
                    ss << std::setw( 2 ) << std::setfill( '0' ) << hours << ":";
                    ss << std::setw( 2 ) << std::setfill( '0' ) << minutes << ":";
                    ss << std::setw( 2 ) << std::setfill( '0' ) << seconds;
                
                    output += prefixSeparator + ss.str() + suffixSeparator;
                }
            
            }
        
            RBOUT(output);
        
        }
        
    } // end if enabled
    
    currentGen = gen;
}


/** Print header for monitored values */
void ScreenMonitor::printHeader( void )
{
    // only print header for cold chains and the enable replicate
    if ( enabled == true && mcmc->getChainPosteriorHeat() == 1.0)
    {
        // print empty line first
        RBOUT("\n");
    
        // print everything to a string stream
        std::stringstream ss;

        // print one column for the iteration number
        std::string header = "Iter";

        int width = 9;
    
        int numWidth = int( log10( numCycles ) ) + 1;
        width = width > numWidth ? width : numWidth;
    
        int columnWidth = 12;

        StringUtilities::fillWithSpaces( header, width, true );
        ss << header << suffixSeparator;
    
        if ( printChainIndex )
        {
            header = "Chain Index";
            StringUtilities::fillWithSpaces( header, columnWidth, false );
            ss << prefixSeparator << header << suffixSeparator;
        }
        
        if ( posterior )
        {
            header = "Posterior";
            StringUtilities::fillWithSpaces( header, columnWidth, false );
            ss << prefixSeparator << header << suffixSeparator;
        }
    
        if ( likelihood )
        {
            header = "Likelihood";
            StringUtilities::fillWithSpaces( header, columnWidth, false );
            ss << prefixSeparator << header << suffixSeparator;
        }
    
        if ( prior )
        {
            header = "Prior";
            StringUtilities::fillWithSpaces( header, columnWidth, false );
            ss << prefixSeparator << header << suffixSeparator;
        }
    
        for (std::vector<DagNode *>::const_iterator it=nodes.begin(); it!=nodes.end(); it++)
        {
            const DagNode* theNode = *it;

            // print the header
            if ( theNode->getName() != "" )
            {
                ss << prefixSeparator;
                theNode->printName( ss, prefixSeparator + suffixSeparator, int( columnWidth ), false );
                ss << suffixSeparator;
            }
            else
            {
                header = "Unnamed";
                StringUtilities::fillWithSpaces( header, columnWidth, false );
                ss << prefixSeparator << header << suffixSeparator;
            }
        }
        
        if ( printElapsedTime )
        {
            // We know it takes 8 characters to print the waiting time, so hard-set column
            // width to this number
            header = "elapsed";
            StringUtilities::fillWithSpaces( header, 8, false );
            ss << prefixSeparator << header << suffixSeparator;
        }
    
        if ( printWaitingTime )
        {
            // We know it takes 8 characters to print the waiting time, so hard-set column
            // width to this number
            header = "ETA";
            StringUtilities::fillWithSpaces( header, 8, false );
            ss << prefixSeparator << header << suffixSeparator;
        }

        RBOUT(ss.str());
        
        std::string dashes = "";
        for (size_t i=0; i<ss.str().size(); ++i)
        {
            dashes += "-";
        }
        RBOUT(dashes);
        
    } // end if enabled
    
}


void ScreenMonitor::reset( size_t n )
{
    startGen = currentGen;
    numCycles = n;
    startTime = time( NULL );
    printWaitingTime = true;
}



/**
 * Set the replicate index.
 * If the index is not 1, then the monitor will be disabled.
 */
void ScreenMonitor::setReplicateIndex(size_t idx)
{
    
    // store the index for possible later uses
    replicateIndex = idx;
    
    enabled = replicateIndex == 1;
    
}


