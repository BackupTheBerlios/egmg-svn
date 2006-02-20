/** \file Prolongation.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface of the abstract class Prolongation
 */
#ifndef PROLONGATION_H_
#define PROLONGATION_H_
#include<valarray>
#include "../general/parameters.h"
#include "../general/Transfer_Operator.h"

namespace mg
{
	
//forward declaration because Stencil uses Prolongation
class Stencil;

/**
 * \brief Prolongation is the interface of a 2D Prolongation operator
 */
class Prolongation : public Transfer_Operator
{
public:
	virtual ~Prolongation() {};
	/**
	 * \brief prolong does an interpolation on the input vector
	 * 
	 * \param[in] u			the vector that represents a rectangel to prolongate
	 * \param[in] stencil	the stencil rep. of the pde needed for matrix dep.
	 * 						Prolongations.
	 * \param[in] Nx		Number of steps in x direction
	 * \param[in] Ny		Number of steps in y direction
	 * \return 				a vector representing the prolongated rectangel
	 */
	virtual std::valarray<precision> prolong(const std::valarray<precision>& u,
									const Stencil& stencil,
									const size_t Nx,const size_t Ny) const =0;
};

}

#endif /*PROLONGATION_H_*/
