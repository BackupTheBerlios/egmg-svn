/** \file Biharmonic2D2.h
 * \author André Oeckerath
 * \brief Contains the class Biharmonic2D2
 */
#ifndef BIHARMONIC2D2_H_
#define BIHARMONIC2D2_H_

#include <vector>
#include "Stencil.h"

namespace mg
{

/**
 * \brief   Biharmonic2D2 represents a discrete laplace like operator of fourth
 *          order.
 * 
 * Biharmonic2D2 is the stencil representing the discrete differential operator
 * \f[
 *  L_hu_h := (\delta_{xx} + \delta_{yy})(\delta_{xx} + \delta_{yy}) u_h
 * \f]
 */
class Biharmonic2D2 : public Stencil
{
public:
    /**
     * \brief The constructor of a Biharmonic2D2 object
     */
    explicit Biharmonic2D2();
    virtual ~Biharmonic2D2();

    Precision apply(
        const NumericArray& u,
        const Position pos,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny ) const;

    Precision getCenter(
        const Position pos,
        const Index,
        const Index,
        const Index nx,
        const Index ny ) const;

    /**
     * \brief retruns the coefficents of Biharmonic2D2 for a center point.
     * 
     * get_L_c() returns the coefficients of Biharmonic2D2 as valarray. For the
     * ordering of the elements
     * \see Stencil
     * E.g. Biharmonic2D2 with the stepsize hx=hy=1 at a center point
     * looks like:
     * \f[
     * \frac{1}{12}
     * \left[\begin{array}{ccccc}
     *  &   &       1   &   &   \\
     *  &   2&      -8  2&  &   \\
     *  1&  -8& 20& -8& 1\\
     *  &   2&      -8  2&  &   \\
     * &    &       1   &   &   \\
     * \end{array}\right]
     * \f]
     * So we have:\n
     * L   ={20,-8,-8,-8,-8, 1, 1, 1, 1, 2, 2, 2, 2}\n
     * J_x ={ 0,-1,0,1,0,-2,0,2,0,-1,1,1,-1}\n
     * J_y ={ 0,0,1,0,-1,0,2,0,-2,1,1,-1,-1}\n
     * 
     * \param[in] pos   relative position in the domain according to the
     *                  enum Postion 
     * \param[in]       the x coordinate of the center element (not used)
     * \param[in]       the y coordinate of the center element (not used)
     * \param[in] nx    the step size in x direction
     * \param[in] ny    the step size in y direction
     * \return          the coefficients of Biharmonic2D4
     */
    NumericArray getL(
        const Position pos,
        const Index,
        const Index,
        const Index nx,
        const Index ny ) const;

	PositionArray getJx(
		const Position pos,
		const Index,
		const Index ) const;

    PositionArray getJy(
		const Position pos,
		const Index,
		const Index ) const;
    
    /**
     * \brief gives the max expansion of Biharmonic2D2
     * 
     * \return  2
     */
    Index size() const;

    /**
     * \brief returns true, because Biharmonic2D2 is constant
     * 
     * \return  true
     */
    bool isConstant() const;

private:
	std::vector<PositionArray > initJx_();
	std::vector<PositionArray > initJy_();
	std::vector<PositionArray > jx_;
	std::vector<PositionArray > jy_;

};

}

#endif /*BIHARMONIC2D2_H_*/
