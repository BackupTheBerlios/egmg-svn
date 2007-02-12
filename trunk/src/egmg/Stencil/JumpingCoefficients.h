/** \file JumpingCoefficients.h
 * \author Kai Ruhnau, Jiri Kraus, Roman Wienands
 * \brief Contains the class JumpingCoefficients
 */
#ifndef JUMPINGCOEFFICIENTS_H_
#define JUMPINGCOEFFICIENTS_H_

#include "Stencil.h"

namespace mg
{

/**
 * \brief   A stencil with jumping coefficients.
 * 
 * JumpingCoefficients is the stencil representing the discrete differential operator
 * \f[
 * follows
 * \f]
 */
class JumpingCoefficients : public Stencil
{
private:
    const PositionArray jx_;
    const PositionArray jy_;
    int mode_;
    
    //initilize jx_, makes it possible to make jx_ const
    PositionArray initJx_() const;
    //initilize jy_, makes it possible to make jy_ const
    PositionArray initJy_() const;
    
    Precision a_(
        const Precision x,
        const Precision y) const;
    
    Precision b_(
        const Precision x,
        const Precision y) const;
public:
	enum Mode { mA, mB, mC, mD };
    /**
     * \brief The constructor of a JumpingCoefficients object
     * 
     * JumpingCoefficients constructs a JumpingCoefficients
     * 
     * \param[in] mode      (default mA)
     */
    explicit JumpingCoefficients(Mode mode = mA );

    virtual ~JumpingCoefficients();

    Precision apply(
        const NumericArray& u,
        const Position,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny) const;

    Precision getCenter(
       const Position,
       const Index sx,
       const Index sy,
       const Index nx,
       const Index ny) const;

    NumericArray getL(
        const Position,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny) const;

    PositionArray getJx(
		const Position,
		const Index,
		const Index ) const;

    PositionArray getJy(
		const Position,
		const Index,
		const Index ) const;

    /**
     * \brief gives the max expansion of JumpingCoefficients
     * 
     * \return  1
     */
	Index size() const;

    /**
     * \brief returns false, because JumpingCoefficients is not constant
     * 
     * \return  false
     */
    bool isConstant() const;
};

}

#endif /*JUMPINGCOEFFICIENTS_H_*/
