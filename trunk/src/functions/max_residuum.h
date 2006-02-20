/** \file max_residuum.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface of the function max_residuum
 */
#ifndef MAX_RESIDUUM_H_
#define MAX_RESIDUUM_H_
#include "../general/parameters.h"
#include "../Stencil/Stencil.h"

namespace mg
{
	
	/**
	 * \brief max_residuum calculates the maximum norm of the residuum
	 * 
	 * max_residuum calcutates the maximum norm of the residuum (defect) of the
	 * discrete PDE given by stencil, i.e
	 * \f[
	 * \parallel d_h\parallel_\infty =\parallel f_h-L_h\cdot u_h\parallel_\infty
	 * \f]
	 * 
	 * \param[in] u			the acctual approximation of the solution
	 * \param[in] fv		the right hand side of the pde
	 * \param[in] stencil	the stencil rep. of the pde
	 * \param[in] Nx		number of steps in x direction
	 * \param[in] Ny		number of steps in y direction
	 * \return				the maximum norm of the residuum
	 */
	precision max_residuum(const std::valarray<precision>& u,
						const std::valarray<precision>& fv,
						const Stencil& stencil,
						const size_t Nx,const size_t Ny);
}

#endif /*MAX_RESIDUUM_H_*/
