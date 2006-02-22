/** \file twonorm_residuum.h
 * \author Matthias Rettenmeier
 * \brief contains the interface of the function twonorm_residuum
 */
#ifndef TWONORM_RESIDUUM_H_
#define TWONORM_RESIDUUM_H_

#include <valarray>
#include "../general/parameters.h"
#include "../Stencil/Stencil.h"

namespace mg
{
    
/**
 * \brief twonorm_residuum() calculates the two norm of the residuum
 * 
 * twonorm_resid() is a memberfunction of the relaxation interface to avoid
 * using two diffrent discretizations by accident.
 * 
 * \param u     the acctual approximation of the solution
 * \param fv        the right hand side of the pde
 * \param nx    number of steps in x direction
 * \param ny    number of steps in y direction
 * \return      the two norm of the resduum
 */
    Precision twonorm_residuum(
        const std::valarray<Precision>& u,
        const std::valarray<Precision>& fv,
        const Stencil& stencil,
        const size_t nx,
        const size_t ny);
}

#endif /*TWONORM_RESIDUUM_H_*/
