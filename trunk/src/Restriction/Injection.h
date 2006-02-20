/** \file Injection.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface of the class Injection
 */
#ifndef INJECTION_H_
#define INJECTION_H_

#include "Restriction.h"



 
namespace mg
{

/**
 * \brief Injection is a 2D restriction operator
 * 
 * Injection is a representation of a 2D restriction operator that uses
 * injection to do its job.
 */
class Injection : public mg::Restriction
{
private:
	const precision weight_;
	const std::valarray<precision> I;
	const std::valarray<int> J_x;
	const std::valarray<int> J_y;
	//initilize I, makes it possible to make I const
	std::valarray<precision> init_I(precision weight) const
	{
		return std::valarray<precision>(weight,1);
	}
	//initilize J_x, makes it possible to make J_x const
	std::valarray<int> init_J_x() const
	{
		return std::valarray<int>(0,1);
	}
	//initilize J_y, makes it possible to make J_y const
	std::valarray<int> init_J_y() const
	{
		return std::valarray<int>(0,1);
	}
	//we don't want these autogenerated contors and operators
	Injection(const Injection& rhs);
	Injection& operator=(const Injection& rhs);
public:
	/**
	 * \brief The constructor of a Injection object
	 * 
	 * Injection constructs a Injection object with:
	 * \param[in] weight	the weight to do injection with (default 1.0)
	 */
	Injection(const precision weight =1.0)
		: weight_(weight),I(init_I(weight)), J_x(init_J_x()), J_y(init_J_y()) {}
	virtual ~Injection(){}
	/**
	 * \brief restriction() restricts the given vector to a smaler gird
	 * 
	 * restriction() restricts the given vector which represents a rectangular
	 * grid to a smaler grid whith the double step size with weighted injection.
	 * The weight is set with the constructor
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
	std::valarray<precision> restriction(
					const std::valarray<precision>& u, const Stencil& stencil,
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

#endif /*INJECTION_H_*/
