/** \file directLUSolver.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief directLUSolver.h contains the declaration of the func. directLUSolver.
 */
#ifndef DIRECTLUSOLVER_H_
#define DIRECTLUSOLVER_H_

#include "../general/parameters.h"
#include "../Stencil/Stencil.h"

namespace mg
{

/**
 * \brief directLUSolver() solves the discrete pde given by stencil
 * 
 * directLUSolver() solves the discrete linear pde given by stencil with
 * LU decomposition.
 * 
 * \param[in,out] u     the solution of the discrete pde
 * \param[in] f         the right hand side of the discrete pde
 * \param[in] stencil   the stencil rep. of the discrete pde
 * \param[in] nx        number of steps in x direction
 * \param[in] ny        number of steps in y direction
 */
void directLUSolver(
    NumericArray& u,
    const NumericArray& f,
    const Stencil& stencil,
    const Index nx,
    const Index ny);
}
#endif /*DIRECTLUSOLVER_H_*/
