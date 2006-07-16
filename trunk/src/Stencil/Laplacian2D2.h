/** \file Laplacian2D2.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the class Laplacian2D2
 */
#ifndef LAPLACIAN2D2_H_
#define LAPLACIAN2D2_H_

#include "Stencil.h"

namespace mg
{

/**
 * \brief   Laplacian2D2 represents a discrete laplace like operator of second
 *          order.
 * 
 * Laplacian2D2 is the stencil representing the discrete differential operator
 * \f[
 *  L_hu_h := a_x(u_h)_{xx}+a_y(u_h)_{yy}
 * \f]
 */
class Laplacian2D2 : public Stencil
{
public:
    /**
     * \brief The constructor of a Laplacian2D2 object
     * 
     * Laplacian2D2 constructs a Laplacian2D2 where \f$a_x\f$ and \f$a_y\f$
     * are given by:
     * \param[in] ax   coefficient of the diff. operator (default 1.0)
     * \param[in] ay   coefficient of the diff. operator (default 1.0)
     */
    explicit Laplacian2D2(
        Precision ax = 1.0,
        Precision ay = 1.0 );

    virtual ~Laplacian2D2();

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
     * \brief gives the max expansion of Laplacian2D2
     * 
     * \return  1
     */
    Index size() const;

    /**
     * \brief returns true, because Laplacian2D2 is constant
     * 
     * \return  true
     */
    bool isConstant() const;
private:
    const PositionArray jx_;
    const PositionArray jy_;
    const Precision ax_;
    const Precision ay_;
    //initilize jx_, makes it possible to make jx_ const
    PositionArray initJx_() const;
    //initilize jy_, makes it possible to make jy_ const
    PositionArray initJy_() const;
};

}

#endif /*LAPLACIAN2D2_H_*/
