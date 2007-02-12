/** \file MSV2D4.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the class MSV2D4
 */
#ifndef MSV2D4_H_
#define MSV2D4_H_

#include "Stencil.h"
#include "../Prolongation/Prolongation.h"
#include "../Restriction/Restriction.h"

namespace mg
{

/**
 * \brief   MSV2D4 represents a discrete laplace like operator of second
 *          order.
 * 
 * MSV2D4 is the stencil representing the discrete differential operator
 * \f[
 *  L_h u_h := a_x(u_h)_{xx}+a_y(u_h)_{yy}
 * \f]
 */
class MSV2D4 : public Stencil
{
public:
    /**
     * \brief The constructor of a MSV2D4 object
     * 
     * MSV2D4 constructs a MSV2D4 where \f$a_x\f$ and \f$a_y\f$
     * are given by:
     * \param[in] ax   coefficient of the diff. operator (default 1.0)
     * \param[in] ay   coefficient of the diff. operator (default 1.0)
     */
    explicit MSV2D4( Precision ax = 1.0, Precision ay = 1.0);

    virtual ~MSV2D4();

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
     * \brief gives the max expansion of MSV2D4
     * 
     * \return  1
     */
    Index size() const;

    /**
     * \brief returns true, because MSV2D4 is constant
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
