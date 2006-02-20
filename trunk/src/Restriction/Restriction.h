/** \file Restriction.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface Restriction
 */
#ifndef RESTRICTION_H_
#define RESTRICTION_H_
#include<valarray>
#include "../general/parameters.h"
#include "../Prolongation/Prolongation.h"

namespace mg
{
	
//forward declaration because Stencil uses Restriction
class Stencil;

/**
 * \brief Restriction is a 2D restriction operator
 */
class Restriction : public mg::Transfer_Operator
{
public:
	virtual ~Restriction() {}
	/**
	 * \brief restriction() restricts the given vector to a smaller grid
	 * 
	 * restriction() restricts the given vector which represents a rectangular
	 * grid to a smaller grid.
	 * \param[in] u					the vector to restrict
	 * \param[in] stencil			the stencil rep. of the pde needed
	 * 								for matrix dep. Restrictions
	 * \param[in] prolong			Prolongation used, needed for matrix dep.
	 * 								Restrictions
	 * \param[in] Nx				number of steps in x direction
	 * \param[in] Ny				number of steps in y direction
	 * \throw std::domain_error		if Nx or Ny is not divedable by 2
	 * \return						a vector with the values on the restricted
	 * 								grid
	 */
	virtual std::valarray<precision> restriction(
					const std::valarray<precision>& u, const Stencil& stencil,
					const Prolongation& prolong, 
					const size_t Nx, const size_t Ny) const =0;
};

}

#endif /*RESTRICTION_H_*/
