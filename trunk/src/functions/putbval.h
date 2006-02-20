/** \file putbval.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface of the function putbval
 */
#ifndef PUTBVAL_H_
#define PUTBVAL_H_
#include "../general/parameters.h"

namespace mg
{
	/**
	 * \brief putbval() fills the given vector with diriclet boundry conditions
	 * 
	 * \param[in,out] u		the vector to be filled on the boundries of the
	 * 						rectangular grid with g
	 * \param[in] Nx		number of steps in x direction
	 * \param[in] Ny		number of steps in y direction
	 * \param[in] g			the diriclet boundry condition
	 */
	void putbval(std::valarray<precision>& u,const size_t Nx, const size_t Ny,
					const function2D g);
}

#endif /*PUTBVAL_H_*/
