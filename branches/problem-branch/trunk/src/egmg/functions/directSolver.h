/** \file directSolver.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 *  <a href="mailto:kai.lists.berlios.de@tragetaschen.dyndns.org">Kai Ruhnau</a>
 * \brief directSolver.h contains the declaration of the func. directSolver.
 */
#ifndef DIRECTSOLVER_H_
#define DIRECTSOLVER_H_

#include "../Problem/Problem.h"

namespace mg
{

/**
 * \brief directSolver() solves the discrete pde given by stencil
 * 
 * directSolver() solves the discrete linear pde given by stencil with
 * LU decomposition.
 * 
 * \param[in,out] u     the solution of the discrete pde
 * \param[in] f         the right hand side of the discrete pde
 * \param[in] stencil   the stencil rep. of the discrete pde
 * \param[in] nx        number of steps in x direction
 * \param[in] ny        number of steps in y direction
 */
NumericArray directSolver(
    const LinearEquationSystem& les);
}
#endif /*DIRECTSOLVER_H_*/
