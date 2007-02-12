/** \file Helmholtz2D2.h
 * \author Andre Oeckerath
 * \brief Contains the class Helmholtz2D2
 */
#ifndef Helmholtz2D2_H_
#define Helmholtz2D2_H_

#include "Stencil.h"

namespace mg
{

/**
 * \brief   Helmholtz2D2 represents a discrete laplace like operator of second
 *          order.
 * 
 * Helmholtz2D2 is the stencil representing the discrete differential operator
 * \f[
 *  L_h u_h := a_x(u_h)_{xx}+a_y(u_h)_{yy}+c_{xy}(u_h)
 * \f]
 */
class Helmholtz2D2 : public Stencil
{

public:
    /**
     * \brief The constructor of a Helmholtz2D2 object
     * 
     * Helmholtz2D2 constructs a Helmholtz2D2 where \f$a_x\f$ , \f$a_y\f$
     * and \f$c_{xy}\f$ are given by:
     * \param[in] ax   coefficient of the diff. operator (default 1.0)
     * \param[in] ay   coefficient of the diff. operator (default 1.0)
     * \param[in] c coefficient of the diff. operator (default 1.0)
     */
    explicit Helmholtz2D2(
        Precision ax = 1.0,
        Precision ay = 1.0,
        Precision c = 1.0 );

    virtual ~Helmholtz2D2();

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
     * \brief gives the max expansion of Helmholtz2D2
     * 
     * \return  1
     */
    Index size() const;

    /**
     * \brief returns true, because Helmholtz2D2 is constant
     * 
     * \return  true
     */
    bool isConstant() const;
private:
    const PositionArray jx_;
    const PositionArray jy_;
    const Precision ax_;
    const Precision ay_;
    const Precision c_;
    //initilize jx_, makes it possible to make jx_ const
    PositionArray initJx_() const;
    //initilize jy_, makes it possible to make jy_ const
    PositionArray initJy_() const;

};

}

#endif /*Helmholtz2D2_H_*/
