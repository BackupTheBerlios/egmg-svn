/** \file twonormResiduum.h
 * \author Matthias Rettenmeier
 * \brief contains the interface of the function twonormResiduum
 */
#ifndef TWONORMRESIDUUM_H_
#define TWONORMRESIDUUM_H_


#include "../general/parameters.h"
#include "../Stencil/Stencil.h"

namespace mg
{
    
/**
 * \brief twonormResiduum() calculates the two norm of the residuum
 *
 * \param[in] u         the acctual approximation of the solution
 * \param[in] f         the right hand side of the pde
 * \param[in] stencil   the stencil rep. of the pde
 * \param[in] nx        number of steps in x direction
 * \param[in] ny        number of steps in y direction
 * \return              the two norm of the resduum
 */
    Precision twonormResiduum(
        const NumericArray& u,
        const NumericArray& f,
        const Stencil& stencil,
        const Index nx,
        const Index ny);
}

#endif /*TWONORMRESIDUUM_H_*/
