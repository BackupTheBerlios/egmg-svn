/** \file twonorm_residuum.h
 * \author Matthias Rettenmeier
 * \brief contains the interface of the function twonorm_residuum
 */
#ifndef TWONORM_RESIDUUM_H_
#define TWONORM_RESIDUUM_H_
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
 * \param u		the acctual approximation of the solution
 * \param fv		the right hand side of the pde
 * \param Nx	number of steps in x direction
 * \param Ny	number of steps in y direction
 * \return		the two norm of the resduum
 */
precision twonorm_residuum(const std::valarray<precision>& u,
						const std::valarray<precision>& fv,
						const Stencil& stencil,
						const size_t Nx,const size_t Ny);
}

#endif /*TWONORM_RESIDUUM_H_*/
