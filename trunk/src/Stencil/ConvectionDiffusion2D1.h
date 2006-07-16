/** \file ConvectionDiffusion2D1.h
 * \author Andre Oeckerath
 * \brief Contains the class ConvectionDiffusion2D1
 */
#ifndef CONVECTIONDIFFUSION2D1_H_
#define CONVECTIONDIFFUSION2D1_H_

#include "Stencil.h"

namespace mg
{

/**
 * \brief   ConvectionDiffusion2D1 represents a discrete Convection Diffusion operator of first
 *          order.
 * 
 * ConvectionDiffusion2D1 is the stencil representing the discrete differential operator
 * \f[
 *  L_h u_h := \epsilon (-(u_h)_{xx}-(u_h)_{yy}) + a(x/y) u_x + b(x,y) u_y
 * \f]
 */
class ConvectionDiffusion2D1 : public Stencil
{
private:
    const PositionArray jx_;
    const PositionArray jy_;
    const Precision epsilon_;
    int mode_;
    
    //initilize jx_, makes it possible to make jx_ const
    PositionArray initJx_() const;
    //initilize jy_, makes it possible to make jy_ const
    PositionArray initJy_() const;
    
    Precision a_(
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny) const;
    
    Precision b_(
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny) const;
public:
	enum Mode { A, B };
    /**
     * \brief The constructor of a ConvectionDiffusion2D1 object
     * 
     * ConvectionDiffusion2D1 constructs a ConvectionDiffusion2D1 where
     * \f$\epsilon\f$ is given by:
     * 
     * \param[in] epsilon   coefficient of the diffusion part (default 1)
     * \param[in] mode      (default A)
     */
    explicit ConvectionDiffusion2D1(Precision epsilon = 1.0, Mode mode = A );

    virtual ~ConvectionDiffusion2D1();

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
     * \brief gives the max expansion of ConvectionDiffusion2D1
     * 
     * \return  1
     */
	Index size() const;

    /**
     * \brief returns false, because ConvectionDiffusion2D1 is not constant
     * 
     * \return  false
     */
    bool isConstant() const;
};

}

#endif /*CONVECTIONDIFFUSION2D1_H_*/
