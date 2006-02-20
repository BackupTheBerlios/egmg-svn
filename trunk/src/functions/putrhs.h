/** \file putrhs.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface of the function putrhs
 */
#ifndef PUTRHS_H_
#define PUTRHS_H_
#include "../general/parameters.h"

namespace mg
{
	/**
	 * \brief putrhs() fills the given vector with the rhs of the pde
	 * 
	 * \param[in,out] fv	the vector of the right hand side of the pde
	 * \param[in] Nx		number of steps in x direction
	 * \param[in] Ny		number of steps in y direction
	 * \param[in] f			the right hand side function of the pde	
	 */
	void putrhs(std::valarray<precision>& fv,const size_t Nx, const size_t Ny,
					const function2D f);
}

#endif /*PUTRHS_H_*/
