/** \file ConDiff2D1.h
 * \author Andre Oeckerath
 * \brief Contains the class ConDiff2D1
 */
#ifndef CONDIFF2D1_H_
#define CONDIFF2D1_H_

#include "Stencil.h"

namespace mg
{

/**
 * \brief   ConDiff2D1 represents a discrete Convection Diffusion operator of
 *          first order.
 * 
 * ConDiff2D1 is the stencil representing the discrete differential operator
 * \f[
 *  L_h u_h := \epsilon (-(u_h)_{xx}-(u_h)_{yy})
 *             + cos(\beta) u_x + sin(\beta) u_y
 * \f]
 */
class ConDiff2D1 : public Stencil
{
private:
    const PositionArray jx_;
    const PositionArray jy_;
    const Precision epsilon_;
    const Precision a1_;
    const Precision a2_;
    //initilize jx_, makes it possible to make jx_ const
    PositionArray initJx_() const;
    //initilize jy_, makes it possible to make jy_ const
    PositionArray initJy_() const;

public:
    /**
     * \brief The constructor of a ConDiff2D1 object
     * 
     * ConDiff2D1 constructs a ConDiff2D1 where \f$\epsilon\f$, \f$a_1\f$ and
     * \f$a_2\f$ are given by:\n
     * \f[
     * a_1 = cos(\beta)
     * \f]
     * \f[
     * a_2 = sin(\beta)
     * \f]
     * \param[in] epsilon coefficient of the diffusion part (default 1.0)
     * \param[in] beta    coefficient of the convection part (default pi/4)
     */
    explicit ConDiff2D1(
        Precision epsilon = 1.0,
        Precision beta = 0.78539816339745 );

    virtual ~ConDiff2D1();
    
    Precision apply(
        const NumericArray& u,
        const Position,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny) const;

    Precision getCenter(
        const Position,
        const Index,
        const Index,
        const Index nx,
        const Index ny) const;

    NumericArray getL(
        const Position,
        const Index,
        const Index,
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
     * \brief gives the max expansion of ConDiff2D1
     * 
     * \return  1
     */
    Index size() const;

    /**
     * \brief returns true, because ConDiff2D1 is constant
     * 
     * \return  true
     */
    bool isConstant() const;
    
};

}

#endif /*CONDIFF2D1_H_*/
