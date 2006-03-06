/** \file output.h
 * \author Matthias Rettenmeier
 * \brief contains the interface of the functions onvergenceRates
 * \todo implement a gnuplot output function
 */
#ifndef OUTPUT_H_
#define OUTPUT_H_


#include<vector>
#include "../general/parameters.h"

namespace mg
{

/**
 * \brief convergenceRates prints convergence rates to stdout
 * 
 * Upon the a history of residues convergenceRates calculatetes convergence 
 * rates and outputs them.
 * CRI = ( ||x{i}||/||x{0}||)^(1/i)
 * CR5 = ( ||x{i}||/||x{5}||)^(1/(i-5)) for i>5
 * CR  = ( ||x{i}||/||x{i-1}||
 * \param[in] vec    a vector containing a history of residues messured 
 *                   in the same norm
 */
void convergenceRates(const std::vector<Precision>& vec, std::ostream& out);

/**
 * 
 */
void gnuPlotDiscreteFunction(
    const NumericArray& u,
    const Index nx,
    const Index ny,
    std::ostream& out);

}

#endif /*OUTPUT_H_*/
