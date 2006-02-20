/** \file residuum.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface of the function residuum
 */
#ifndef RESIDUUM_H_
#define RESIDUUM_H_
#include "../general/parameters.h"
#include "../Stencil/Stencil.h"

namespace mg
{
	/**
	 * \brief residuum calculates the residuum
	 * 
	 * residuum calcutates the residuum (defect) of the discrete PDE given
	 * by stencil, i.e
	 * \f[ d_h = f_h - L_h\cdot u_h \f]
	 * 
	 * \param[in] u			the acctual approximation of the solution
	 * \param[in] fv		the right hand side of the pde
	 * \param[in] stencil	the stencil rep. of the pde
	 * \param[in] Nx		number of steps in x direction
	 * \param[in] Ny		number of steps in y direction
	 * \return				a vector with the resduum
	 */
	std::valarray<precision> residuum(const std::valarray<precision>& u,
								   const std::valarray<precision>& fv,
								   const Stencil& stencil,
								   const size_t Nx,const size_t Ny);
}

#endif /*RESIDUUM_H_*/
