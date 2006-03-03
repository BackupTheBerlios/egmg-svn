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
 * \param u     the acctual approximation of the solution
 * \param f     the right hand side of the pde
 * \param nx    number of steps in x direction
 * \param ny    number of steps in y direction
 * \return      the two norm of the resduum
 */
    Precision twonormResiduum(
        const NumericArray& u,
        const NumericArray& f,
        const Stencil& stencil,
        const size_t nx,
        const size_t ny);
}

#endif /*TWONORMRESIDUUM_H_*/
