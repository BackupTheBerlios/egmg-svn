/** \file Half_Weighting.h
 * \author Mareike Riepl, Michaela Anghel, Daniela Steffes-Lai
 * \brief Contains the interface of the class Half Weighting
 */
#ifndef HALF_WEIGHTING_H_
#define HALF_WEIGHTING_H_

#include "Restriction.h"

namespace mg
{

/**
 * \brief Half_Weighting is a 2D restriction operator
 * 
 * Half_Weighting is a representation of a 2D restriction operator that uses
 * half weighting to do its job.
 */
class Half_Weighting : public mg::Restriction
{
private:
	const precision weight_;
	const std::valarray<precision> I;
	const std::valarray<int> J_x;
	const std::valarray<int> J_y;
	//initilize I, makes it possible to make I const
	std::valarray<precision> init_I(precision weight) const
	{
		const precision t[] = {weight/4,weight/8,weight/8,weight/8,weight/8,
									0,0,0,0};
		return std::valarray<precision>(t,9);
	}
	//initilize J_x, makes it possible to make J_x const
	std::valarray<int> init_J_x() const
	{
		const int t[] = {0,-1,0,1,0,-1,1,1,-1};
		return std::valarray<int>(t,9);
	}
	//initilize J_y, makes it possible to make J_y const
	std::valarray<int> init_J_y() const
	{
		const int t[] = {0,0,1,0,-1,1,1,-1,-1};
		return std::valarray<int>(t,9);
	}
	//we don't want these autogenerated contors and operators
	Half_Weighting(const Half_Weighting& rhs);
	Half_Weighting& operator=(const Half_Weighting& rhs);
public:
	/**
	 * \brief The constructor of a Half_Weighting object
	 * 
	 * Half_Weighting constructs a Half_Weighting object with:
	 * \param[in] weight	the weight to do half weighting with (default 1.0)
	 */
	Half_Weighting(precision weight = 1.0)
		: weight_(weight),I(init_I(weight)), J_x(init_J_x()), J_y(init_J_y()) {}
	virtual ~Half_Weighting(){}
	/**
	 * \brief restriction() restricts the given vector to a smaller grid
	 * 
	 * restriction() restricts the given vector which represents a rectangular
	 * grid to a smaller grid with the double step size with the half weighting
	 * Operator.
	 * 
	 * \param[in] u					the vector to restrict
	 * \param[in] stencil			the stencil rep. of the pde needed
	 * 								for matrix dep. Restrictions (not used)
	 * \param[in] prolong			Prolongation used, needed for matrix dep.
	 * 								Restrictions (not used)
	 * \param[in] Nx				number of steps in x direction
	 * \param[in] Ny				number of steps in y direction
	 * \throw std::domain_error		if Nx or Ny is not divedable by 2
	 * \return						a vector with the values on the restricted
	 * 								grid
	 */
	std::valarray<precision> restriction(const std::valarray<precision>& u,
					const Stencil& stencil,
					const Prolongation& prolong,
					const size_t Nx, const size_t Ny) const;
	const std::valarray<precision>& get_I(const size_t, const size_t,
									const size_t, const size_t, const Stencil&) const
	{
		return I;	
	}
	const std::valarray<int>& get_J_x() const
	{
		return J_x;	
	}
	const std::valarray<int>& get_J_y() const
	{
		return J_y;
	}
};

}

#endif /*HALF_WEIGHTING_H_*/
