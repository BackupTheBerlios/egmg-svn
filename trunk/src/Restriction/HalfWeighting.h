/** \file HalfWeighting.h
 * \author Mareike Riepl, Michaela Anghel, Daniela Steffes-Lai
 * \brief Contains the interface of the class HalfWeighting
 */
#ifndef HALFWEIGHTING_H_
#define HALFWEIGHTING_H_

#include "Restriction.h"

namespace mg
{

/**
 * \brief HalfWeighting is a 2D restriction operator
 * 
 * HalfWeighting is a representation of a 2D restriction operator that uses
 * half weighting to do its job.
 */
class HalfWeighting : public mg::Restriction
{
private:
	const Precision weight_;
	const NumericArray i_;
	const PositionArray jx_;
	const PositionArray jy_;
	//initilize i_, makes it possible to make i_ const
	NumericArray initI_(Precision weight) const
	{
		const Precision t[]={weight/4,weight/8,weight/8,weight/8,weight/8,
							 0,0,0,0};
		return NumericArray(t,9);
	}
	//initilize jx_, makes it possible to make jx_ const
	PositionArray initJX_() const
	{
		const int t[]={0,-1,0,1,0,-1,1,1,-1};
		return PositionArray(t,9);
	}
	//initilize jy_, makes it possible to make jy_ const
	PositionArray initJY_() const
	{
		const int t[]={0,0,1,0,-1,1,1,-1,-1};
		return PositionArray(t,9);
	}
	//we don't want these autogenerated ctors and operators
	HalfWeighting(const HalfWeighting& rhs);
	HalfWeighting& operator=(const HalfWeighting& rhs);
public:
	/**
	 * \brief The constructor of a HalfWeighting object
	 * 
	 * HalfWeighting constructs a HalfWeighting object with:
	 * \param[in] weight	the weight to do half weighting with (default 1.0)
	 */
	HalfWeighting(Precision weight=1.0)
		: weight_(weight),
          i_(initI_(weight)),
          jx_(initJX_()), jy_(initJY_()) {}
          
	virtual ~HalfWeighting(){}
    
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
	 * \param[in] prolongation		Prolongation used, needed for matrix dep.
	 * 								Restrictions (not used)
	 * \param[in] nx				number of steps in x direction
	 * \param[in] ny				number of steps in y direction
	 * \throw std::domain_error		if nx or ny is not divedable by 2
	 * \return						a vector with the values on the restricted
	 * 								grid
	 */
	NumericArray restriction(const NumericArray& u,
					const Stencil& stencil,
					const Prolongation& prolongation,
					const Index nx, const Index ny) const;
	const NumericArray& getI(
        const Index,
        const Index,
        const Index,
        const Index,
        const Stencil&) const
	{
		return i_;	
	}
	const PositionArray& getJx() const
	{
		return jx_;	
	}
	const PositionArray& getJy() const
	{
		return jy_;
	}
};

}

#endif /*HALFWEIGHTING_H_*/
