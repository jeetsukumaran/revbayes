#ifndef AbstractRateMatrix_H
#define AbstractRateMatrix_H

#include "MatrixReal.h"
#include "RateMatrix.h"
#include <vector>


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    
    /**
     * @brief Abstract rate matrix class.
     *
     * The abstract rate matrix class provides some basic functionality of most
     * rate matrices. It implements some functionality of the interface RateMatrix.
     * The key element of this abstract class is that it hold a matrix of doubles as the internal value
     * which represent the rates of this matrix. Specific behaviour how the rates are updated,
     * how the stationary frequencies are computed, and how the transition probabilities are computed are
     * left for derived classes.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-07-04, version 1.0
     */
    class AbstractRateMatrix : public RateMatrix {
        
    public:
        virtual                            ~AbstractRateMatrix(void);                                                                   //!< Destructor
        
        // overloaded operators
//        std::vector<double>&                operator[](size_t i);                                                                       //!< Subscript operator
//        const std::vector<double>&          operator[](size_t i) const;                                                                 //!< Subscript operator (const)
        
//        std::vector<std::vector<double> >::const_iterator       begin(void) const;
//        std::vector<std::vector<double> >::iterator             begin(void);
//        std::vector<std::vector<double> >::const_iterator       end(void) const;
//        std::vector<std::vector<double> >::iterator             end(void);
        
        // public methods
        double                              getRate(size_t from, size_t to, double rate=1.0) const;
        double                              getRate(size_t from, size_t to, double age, double rate) const;         //!< Calculate the rate from state i to state j over the given time interval scaled by a rate
        void                                rescaleToAverageRate(double r);                                                             //!< Rescale the rate matrix such that the average rate is "r"
        void                                setDiagonal(void);                                                                          //!< Set the diagonal such that each row sums to zero
        
        // pure virtual methods you have to overwrite
        virtual double                      averageRate(void) const = 0;                                                                //!< Calculate the average rate
        virtual void                        calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const = 0;   //!< Calculate the transition matrix
        virtual AbstractRateMatrix*         clone(void) const = 0;
        virtual const std::vector<double>&  getStationaryFrequencies(void) const = 0;                                                   //!< Return the stationary frequencies
        virtual void                        update(void) = 0;                                                                           //!< Update the rate entries of the matrix (is needed if stationarity freqs or similar have changed)
        

    protected:
        // prevent instantiation
        AbstractRateMatrix(size_t n);                                                                                                   //!< Construct rate matrix with n states
        AbstractRateMatrix(const AbstractRateMatrix& m);                                                                                //!< Copy constructor
        AbstractRateMatrix&                 operator=(const AbstractRateMatrix& r);                                                     //!< Assignment operator
        
        // protected methods available for derived classes
        bool                                checkTimeReversibity(double tolerance);
        
        // protected members available for derived classes
        MatrixReal*                         theRateMatrix;                                                                              //!< Holds the rate matrix
        bool                                needsUpdate;
        
    };
    
}

#endif

