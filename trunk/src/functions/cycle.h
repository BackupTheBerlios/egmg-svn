/** \file cycle.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>, Matthias Rettenmeier
 * \brief Contains the function declaration of cycle
 */
#ifndef CYCLE_H_
#define CYCLE_H_
#include<valarray>
#include "../general/parameters.h"
#include "../Prolongation/Prolongation.h"
#include "../Relaxation/Relaxation.h"
#include "../Restriction/Restriction.h"
#include "../Stencil/Stencil.h"

namespace mg
{
	/**
	 * \brief cycle() does one iterative multigridcycle
	 * 
	 * cycle does one iterative multigridcycle to solve a partial differential
	 * equation (pde) on a discrete rectangular grid. The discretization of the
	 * pde is given by stencil.
	 * Cycle supports V-, F- and W-cycle types. gamma determins what type
	 * of cycle. For gamma = 0 F-cycle is used for gamma>=1 W-cycle with gamma
	 * iterations per grid is used. Trivially V-cycle is used in the case
	 * gamma = 1. Mode is an implementation detail used for recursive calls.
	 * \param[in,out] u				the vector represantation of the unknown
	 * 								function to solve on a rectangular grid.
	 * \param[in] fv				the right hand side of the pde
	 * \param[in] stencil			the stencil representaion of the discrete
	 * 								pde
	 * \param[in] prolong			the prolongation to use
	 * \param[in] restriction		the restriction to use
	 * \param[in] relax				the relaxation to use
	 * \param[in] Nx				number of steps in x direction
	 * \param[in] Ny				number of steps in y direction
	 * \param[in] gamma         	number of iterations to be done at next
	 * 								level (default 1)
	 * \param[in] l             	parameter to set coarsest grid:
	 * 								if min(Nx,Ny)<=2**l coarsest grid reached
	 * 								(default 1)
	 * \param[in] mode          	a flag to remember if f-cycle called
	 * 								(implementaion detail)
	 * \throw std::domain_error		if Nx or Ny is not divedable by 2
	 */
	void cycle(std::valarray<precision>& u,const std::valarray<precision>& f,
                Stencil& stencil,
                const Prolongation& prolong,const Restriction& restriction,
                Relaxation& relax,
                const size_t Nx, const size_t Ny,size_t gamma=1,size_t l=1,
                size_t mode=0);
}

#endif /*CYCLE_H_*/
