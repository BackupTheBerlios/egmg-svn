/** \file output.h
 * \author Matthias Rettenmeier
 * \brief contains the interface of the functions output and outputvalue
 */
#ifndef OUTPUT_H_
#define OUTPUT_H_

#include<valarray>
#include<vector>
#include "../general/parameters.h"

namespace mg
{

/**
 * \brief output generates a basic output
 * Upon the a history of residues output calculatetes convergence rates
 * and outputs them.
 * \param vec    a vector containing a history of residues messured in the same norm
 * \param flag   a flag to let output know which norm is used
 *               1: Max-Norm
 *               2: 2-Norm
 * other norms not yet defined!
 * F1 = ( ||x{i}||/||x{0}||)^(1/i)
 * F2 = ( ||x{i}||/||x{5}||)^(1/(i-5)) for i>5
 * F3 = ( ||x{i}||/||x{i-1}||
 */
int output(std::vector<Precision>& vec, size_t flag);

/**
 * \brief outputvalue outputs a specific point on the grid
 * Upon user input a point on the grid will be outputed if contained on the grid
 * \param u      a valarray containing the calculated solution for the vector u in Au=f.
 * \param Nx	 Number of steps in x direction
 * \param Ny	 Number of steps in y direction
 */
int outputvalue(std::valarray<Precision> u, size_t Nx, size_t Ny);
}

#endif /*OUTPUT_H_*/
