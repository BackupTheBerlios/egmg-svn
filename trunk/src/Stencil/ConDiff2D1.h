/** \file ConDiff2D1.h
 * \author Andr� Oeckerath
 * \brief Contains the class ConDiff2D1
 */
#ifndef CONDIFF2D1_H_
#define CONDIFF2D1_H_

#include <cmath>

#include "Stencil.h"
#include "../Prolongation/Prolongation.h"
#include "../Restriction/Restriction.h"


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
    mutable NumericArray l_;
    const PositionArray jx_;
    const PositionArray jy_;
    const Precision epsilon_;
    const Precision a1_;
    const Precision a2_;
    //initilize jx_, makes it possible to make jx_ const
    PositionArray initJx_() const
    {
        const int t[] = {0,-1,0,1,0};
        return PositionArray(t,5);
    }
    //initilize jy_, makes it possible to make jy_ const
    PositionArray initJy_() const
    {
        const int t[] = {0,0,1,0,-1};
        return PositionArray(t,5);
    }
    //
    //we don't want the autogenerated copy constructor and assignment operator
    ConDiff2D1(const ConDiff2D1&);
    ConDiff2D1& operator=(const ConDiff2D1&);
public:
    /**
     * \brief The constructor of a ConDiff2D1 object
     * 
     * ConDiff2D1 constructs a ConDiff2D1 where \f$\epsilon\f$, \f$a_1\f$ and \f$a_2\f$
     * are given by:
     * \param[in] \epsilon coefficient of the diffusion part (default 1.0)
     * \param[in] a_1 = cos(\beta)  coefficient of the convection part (default \pi/4)
     * \param[in] a_2 = sin(\beta)  coefficient of the convection part (default \pi/4)
     */
    explicit ConDiff2D1(
        Precision epsilon=1.0,
        Precision beta=0.78539816339745) 
            : l_(5),
              jx_(initJx_()), jy_(initJy_()),
              epsilon_(epsilon),
              a1_(std::cos(beta)), a2_(std::sin(beta))
    {}

    virtual ~ConDiff2D1() {}
    inline Precision apply(
        const NumericArray& u,
        const Position,
        const size_t sx,
        const size_t sy,
        const size_t nx,
        const size_t ny) const
    {
        return 
            (2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny
             +fabs(a1_)*nx+fabs(a2_)*ny
            )*u[sy*(nx+1)+sx]
            +(-1.0*epsilon_*nx*nx+(-a1_-fabs(a1_))*nx/2)*u[sy*(nx+1)+sx-1]
            +(-1.0*epsilon_*nx*nx+(a1_-fabs(a1_))*nx/2)*u[sy*(nx+1)+sx+1]
            +(-1.0*epsilon_*ny*ny+(-a2_-fabs(a2_))*ny/2)*u[(sy-1)*(nx+1)+sx]
            +(-1.0*epsilon_*ny*ny+(a2_-fabs(a2_))*ny/2)*u[(sy+1)*(nx+1)+sx];
    }

    inline Precision getCenter(
        const Position,
        const size_t,
        const size_t,
        const size_t nx,
        const size_t ny) const
    {
        return 
            (2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny+fabs(a1_)*nx+fabs(a2_)*ny);
    }

    inline const NumericArray& getL(
        const Position,
        const size_t,
        const size_t,
        const size_t nx,
        const size_t ny) const
    {
        l_[0] = 2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny+fabs(a1_)*nx+fabs(a2_)*ny;
        l_[1] = -1.0*epsilon_*nx*nx+(-a1_-fabs(a1_))*nx/2;
        l_[3] = -1.0*epsilon_*nx*nx+(a1_-fabs(a1_))*nx/2;
        l_[2] = -1.0*epsilon_*ny*ny+(-a2_-fabs(a2_))*ny/2;
        l_[4] = -1.0*epsilon_*ny*ny+(a2_-fabs(a2_))*ny/2;
        return l_;
    }

    inline const PositionArray& getJx(const Position) const
    {
        return jx_;
    }

    inline const PositionArray& getJy(const Position) const
    {
        return jy_;
    }

    /**
     * \brief does nothing for ConDiff2D1
     * \see Stencil
     */
    void pushProlongation(const Prolongation&) {}

    /**
     * \brief does nothing for ConDiff2D1
     * \see Stencil
     */
    void popProlongation() {}

    /**
     * \brief does nothing for ConDiff2D1
     * \see Stencil
     */
    void pushRestriction(const Restriction&) {}

    /**
     * \brief does nothing for ConDiff2D1
     * \see Stencil
     */
    void popRestriction() {}

    /**
     * \brief gives the max expansion of ConDiff2D1
     * 
     * \return  1
     */
    inline size_t size() const
    {
        return 1;
    }

    /**
     * \brief returns true, because ConDiff2D1 is constant
     * 
     * \return  true
     */
    inline bool isConstant() const
    {
        return true;
    }
    
};

}

#endif /*CONDIFF2D1_H_*/
