/** \file Laplacian2D4.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the class Laplacian2D4
 */
#ifndef LAPLACIAN2D4_H_
#define LAPLACIAN2D4_H_

#include <vector>

#include "Stencil.h"

namespace mg
{

/**
 * \brief   Laplacian2D4 represents a discrete laplace like operator of fourth
 *          order.
 * 
 * Laplacian2D4 is the stencil representing the discrete differential operator
 * \f[
 *  L_hu_h := a_x(u_h)_{xx}+a_y(u_h)_{yy}
 * \f]
 */
class Laplacian2D4 : public Stencil
{
public:
    /**
     * \brief The constructor of a Laplacian2D4 object
     * 
     * Laplacian2D4 constructs a Laplacian2D4 where \f$a_x\f$ and \f$a_y\f$
     * are given by:
     * \param[in] ax    coefficient of the diff. operator (default 1.0)
     * \param[in] ay    coefficient of the diff. operator (default 1.0)
     */
    explicit Laplacian2D4( Precision ax = 1.0, Precision ay = 1.0 );

    virtual ~Laplacian2D4();

    Precision apply(
        const NumericArray& u,
        const Position position,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny) const;

    Precision getCenter(
        const Position position,
        const Index,
        const Index,
        const Index nx,
        const Index ny) const;

    /**
     * \brief retruns the coefficents of Laplacian2D4 for a center point.
     * 
     * getL() returns the coefficients of Laplacian2D4 as valarray. For the
     * ordering of the elements
     * \see Stencil
     * E.g. Laplacian2D4 with the stepsize hx=hy=1, a_x=a_y=1.0 at a center point
     * looks like:
     * \f[
     * \frac{1}{12}
     * \left[\begin{array}{ccccc}
     *  &   &       1   &   &   \\
     *  &   &       -16 &   &   \\
     * 1&   -16&    60& -16&    1\\
     *  &   &       -16 &   &   \\
     * &    &       1   &   &   \\
     * \end{array}\right]
     * \f]
     * So we have:\n
     * L    = {60/12,-16/12,-16/12,-16/12,-16/12, 1/12, 1/12,1/12, 1/12}\n
     * J_x  = { 0,-1,0,1,0,-2,0,2,0}\n
     * J_y  = { 0,0,1,0,-1,0,2,0,-2}\n
     * 
     * \param[in]       the x coordinate of the center element (not used)
     * \param[in]       the y coordinate of the center element (not used)
     * \param[in] nx    the step size in x direction
     * \param[in] ny    the step size in y direction
     * \return          the coefficients of Laplacian2D4
     */
    NumericArray getL(
        const Position position,
        const Index,
        const Index,
        const Index nx,
        const Index ny) const;

    PositionArray getJx(
		const Position pos,
		const Index,
		const Index ) const;

    PositionArray getJy(
		const Position pos,
		const Index,
		const Index ) const;

    /**
     * \brief gives the max expansion of Laplacian2D4
     * 
     * \return  2
     */
    Index size() const;

    /**
     * \brief returns true, because Laplacian2D4 is constant
     * 
     * \return  true
     */
    bool isConstant() const;
private:
    const std::vector<PositionArray > jx_;
    const std::vector<PositionArray > jy_;
    const Precision ax_;
    const Precision ay_;

    std::vector<PositionArray > initJx_();
    std::vector<PositionArray > initJy_();
};

}

#endif /*LAPLACIAN2D4_H_*/
