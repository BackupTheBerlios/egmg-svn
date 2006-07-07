/** \file TransposedProlongation.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the interface of the class TransposedProlongation.
 */
#ifndef TRANSPPROLONG_H_
#define TRANSPPROLONG_H_

#include "Restriction.h"
 
namespace mg
{

/**
 * \brief TransposedProlongation is a 2D restriction operator
 * 
 * TransposedProlongation is a representation of a 2D restriction operator that 
 * uses an associated prolongation operator to do its job.
 */
class TransposedProlongation : public mg::Restriction
{
private:
	const Precision weight_;
	mutable NumericArray i_;
	mutable PositionArray jx_;
	mutable PositionArray jy_;
	const Prolongation& prolongation_;

	//we don't want these autogenerated ctors and operators
	TransposedProlongation(const TransposedProlongation& rhs);
	TransposedProlongation& operator=(const TransposedProlongation& rhs);
public:
	/**
	 * \brief The constructor of a TransposedProlongation object
	 * 
	 * TransposedProlongation constructs a TransposedProlongation object with:
	 * \param[in] weight	the weight to do full weighting with (default 1.0)
	 */
	
	TransposedProlongation(
        const Prolongation& prolongation,
        Precision weight = 1.0)
        : weight_(weight), prolongation_(prolongation) {}

	virtual ~TransposedProlongation(){}
    
	/**
	 * \brief restriction() restricts the given vector to a smaller grid
	 * 
	 * restriction() restricts the given vector which represents a rectangular
	 * grid to a smaller grid with the double step size with the associated Prolongation
	 * Operator.
	 * 
	 * \param[in] u					the vector to restrict
	 * \param[in] stencil			the stencil rep. of the pde needed
	 * 								for matrix dep. Restrictions (not used)
	 * \param[in] prolongate		Prolongation used, needed for matrix dep.
	 * 								Restrictions
	 * \param[in] nx				number of steps in x direction
	 * \param[in] ny				number of steps in y direction
	 * \throw std::domain_error		if nx or ny is not divedable by 2
	 * \return						a vector with the values on the restricted
	 * 								grid
	 */
	NumericArray restriction(
        const NumericArray& u,
		const Stencil& stencil,
		const Index nx, const Index ny) const;
        
	const NumericArray& getI(
		const Position pos,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny,
        const Stencil& stencil) const
	{
		Precision scale = 0;
		i_.resize( prolongation_.getI(pos, sx, sy, nx, ny, stencil).size() );
		i_ = prolongation_.getI(pos, sx, sy, nx, ny, stencil);
		for (Index no = 0; no<i_.size(); no++)
			scale += i_[no];
		for (Index no = 0; no<i_.size(); no++)
			i_[no] /= scale;		
		return i_;	
	}
    
	const PositionArray& getJx( const Position pos ) const
	{
		jx_.resize( prolongation_.getJx(pos).size() );
		jx_ = prolongation_.getJx(pos);
		return jx_;	
	}
    
	const PositionArray& getJy( const Position pos ) const
	{
		jy_.resize( prolongation_.getJx(pos).size() );
		jy_ = prolongation_.getJy(pos);
		return jy_;
	}
};

}

#endif /*TRANSPPROLONG_H_*/
