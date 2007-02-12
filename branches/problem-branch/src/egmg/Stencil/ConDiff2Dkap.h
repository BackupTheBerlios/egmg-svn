/** \file ConDiff2Dkap.h
 * \author Andre Oeckerath
 * \brief Contains the class ConDiff2Dkap
 */
#ifndef ConDiff2Dkap_H_
#define ConDiff2Dkap_H_

#include <vector>

#include "Stencil.h"

namespace mg
{

/**
 * \brief   ConDiff2Dkap represents a discrete convection-diffusion operator of 2nd
 *          or 3rd order (depending on kappa).
 * 
 * ConDiff2Dkap is the stencil representing the discrete differential operator
 * \f[
 *  L_h u_h := \epsilon(-(u_h)_{xx}-(u_h)_{yy})+ cos(\beta) u_x + sin(\beta) u_y
 * \f]
 */
class ConDiff2Dkap : public Stencil
{
public:
    /**
     * \brief The constructor of a ConDiff2Dkap object
     * 
     * ConDiff2Dkap constructs a ConDiff2Dkap with:
     * 
     * \param[in] epsilon   (default 1)
     * \param[in] beta      (default pi/4)
     * \param[in] kappa     (default 1/3)
     */
    explicit ConDiff2Dkap(
        Precision epsilon = 1.0,
        Precision beta = 0.78539816339745,
        Precision kappa = 1.0/3 );

    virtual ~ConDiff2Dkap();

    inline Precision apply(
        const NumericArray& u,
        const Position pos,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny) const;

    Precision getCenter(
        const Position pos,
        const Index,
        const Index,
        const Index nx,
        const Index ny) const;

    /**
     * \brief retruns the coefficents of ConDiff2Dkap for a center point.
     * 
     * get_L_c() returns the coefficients of ConDiff2Dkap as valarray. For the
     * ordering of the elements
     * \see Stencil
     * E.g. ConDiv2D4 with the stepsize hx=hy=1, a_1=a_2=1.0 at a center point
     * looks like:
     * e = epsilon_
     * L    = {4e + (1-k)/3, -3(1-k)/4 -e, -(1-k)/4 - e, -(1-k)/4 -e, -3(1-k)/4 -e, (1-k)/4, 0, 0, (1-k)/4}\n
     * jx_  = { 0,-1,0,1,0,-2,0,2,0}\n
     * jy_  = { 0,0,1,0,-1,0,2,0,-2}\n
     * \param[in]       the x coordinate of the center element (not used)
     * \param[in]       the y coordinate of the center element (not used)
     * \param[in] nx    the step size in x direction
     * \param[in] ny    the step size in y direction
     * \return          the coefficients of ConDiff2Dkap
     */
    NumericArray getL(
        const Position pos,
        const Index,
        const Index,
        const Index nx,
        const Index ny) const;

    PositionArray getJx(
		const Position p,
		const Index,
		const Index ) const;

    PositionArray getJy(
		const Position p,
		const Index,
		const Index ) const;

    /**
     * \brief gives the max expansion of ConDiff2Dkap
     * 
     * \return  2
     */
    Index size() const;

    /**
     * \brief returns true, because ConDiff2Dkap is constant
     * 
     * \return  true
     */
    bool isConstant() const;
private:
    const std::vector<PositionArray > jx_;
    const std::vector<PositionArray > jy_;
    const Precision epsilon_;
    const Precision factor_;
    const Precision a1_;
    const Precision a2_;
    
    std::vector<PositionArray > initJx_();
    std::vector<PositionArray > initJy_();
};

}

#endif /*ConDiff2Dkap_H_*/
