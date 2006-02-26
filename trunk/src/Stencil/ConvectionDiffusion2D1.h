/** \file ConvectionDiffusion2D1.h
 * \author Andr� Oeckerath
 * \brief Contains the class ConvectionDiffusion2D1
 */
#ifndef CONVECTIONDIFFUSION2D1_H_
#define CONVECTIONDIFFUSION2D1_H_

#include <cmath>
#include <valarray>
#include <stdexcept>
#include "Stencil.h"
#include "../Prolongation/Prolongation.h"
#include "../Restriction/Restriction.h"

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
    mutable std::valarray<Precision> l_;
    const std::valarray<int> jx_;
    const std::valarray<int> jy_;
    const Precision epsilon_;
    int mode_;
    
    //initilize jx_, makes it possible to make jx_ const
    std::valarray<int> initJx_() const
    {
        const int t[]={0,-1,0,1,0};
        return std::valarray<int>(t,5);
    }
    //initilize jy_, makes it possible to make jy_ const
    std::valarray<int> initJy_() const
    {
        const int t[]={0,0,1,0,-1};
        return std::valarray<int>(t,5);
    }
    
    Precision a_(
        const size_t sx,
        const size_t sy,
        const size_t nx,
        const size_t ny) const
    {
        switch (mode_)
        {
        case 0: return (2.0*sy/ny-1)*(1-(sx*sx)/(nx*nx));
        case 1: return 4.0*sx/nx*(sx/nx-1)*(1-2.0*sy/ny);
        default: ;
        }
    }
    
    Precision b_(
        const size_t sx,
        const size_t sy,
        const size_t nx,
        const size_t ny) const
    {
        switch(mode_)
        {
        case 0: return (2.0*sx*sy/(nx*ny))*(sy/ny-1);
        case 1: return -4.0*sy/ny*(sy/ny-1)*(1-2.0*sx/nx);
        default: ;
        }       
    }
    //we don't want the autogenerated copy constructor  and assignment operator
    ConvectionDiffusion2D1(const ConvectionDiffusion2D1&);
    ConvectionDiffusion2D1& operator=(const ConvectionDiffusion2D1&);
public:
    /**
     * \brief The constructor of a ConvectionDiffusion2D1 object
     * 
     * ConvectionDiffusion2D1 constructs a ConvectionDiffusion2D1 where \f$\epsilon\f$, \f$a_1\f$ and \f$a_2\f$
     * are given by:
     * \param[in] \epsilon coefficient of the diffusion part (default 1.0)
     * \param[in] a=(2 x_2 - 1)(1 - x�_1) coefficient of the convection part (default \pi/4)
     * \param[in] b=2 x_1 x_2 (x_2 - 1)   coefficient of the convection part (default \pi/4)
     */
    explicit ConvectionDiffusion2D1(Precision epsilon =1.0, int mode=0) 
        : l_(5), jx_(initJx_()), jy_(initJy_()), epsilon_(epsilon), mode_(mode)
    {
        if (mode!=0 && mode!=1)
            throw std::range_error("Invalid value for mode");
    }

    virtual ~ConvectionDiffusion2D1() {}

    inline Precision apply(
        const std::valarray<Precision>& u,
        const Position,
        const size_t sx,
        const size_t sy,
        const size_t nx,
        const size_t ny) const
    {
        return 
             (2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny
              +fabs(a_(sx,sy,nx,ny))*nx+fabs(b_(sx,sy,nx,ny))*ny
             )*u[sy*(nx+1)+sx]
            +(-1.0*epsilon_*nx*nx
              +(-a_(sx,sy,nx,ny)-fabs(a_(sx,sy,nx,ny)))*nx/2
             )*u[sy*(nx+1)+sx-1]
            +(-1.0*epsilon_*nx*nx
              +(a_(sx,sy,nx,ny)-fabs(a_(sx,sy,nx,ny)))*nx/2
             )*u[sy*(nx+1)+sx+1]
            +(-1.0*epsilon_*ny*ny
              +(-b_(sx,sy,nx,ny)-fabs(b_(sx,sy,nx,ny)))*ny/2
             )*u[(sy-1)*(nx+1)+sx]
            +(-1.0*epsilon_*ny*ny
              +(b_(sx,sy,nx,ny)-fabs(b_(sx,sy,nx,ny)))*ny/2
             )*u[(sy+1)*(nx+1)+sx];
    }

    inline Precision getCenter(
       const Position,
       const size_t sx,
       const size_t sy,
       const size_t nx,
       const size_t ny) const
    {
        return 
             2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny
            +fabs( a_(sx,sy,nx,ny) )*nx+fabs( b_(sx,sy,nx,ny) )*ny;
    }

    inline const std::valarray<Precision>& getL(
        const Position,
        const size_t sx,
        const size_t sy,
        const size_t nx,
        const size_t ny) const
    {
        l_[0]=
             2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny
            +fabs(a_(sx,sy,nx,ny))*nx+fabs(b_(sx,sy,nx,ny))*ny;
        l_[1]=
             -1.0*epsilon_*nx*nx
            +(-a_(sx,sy,nx,ny)-fabs(a_(sx,sy,nx,ny)))*nx/2;
        l_[3]=
             -1.0*epsilon_*nx*nx
            +(a_(sx,sy,nx,ny)-fabs(a_(sx,sy,nx,ny)))*nx/2;
        l_[2]=
             -1.0*epsilon_*ny*ny
            +(-b_(sx,sy,nx,ny)-fabs(b_(sx,sy,nx,ny)))*ny/2;
        l_[4]=
             -1.0*epsilon_*ny*ny
            +(b_(sx,sy,nx,ny)-fabs(b_(sx,sy,nx,ny)))*ny/2;
        return l_; 
    }

    inline const std::valarray<int>& getJx(const Position) const
    {
        return jx_;
    }

    inline const std::valarray<int>& getJy(const Position) const
    {
        return jy_;
    }

    /**
     * \brief does nothing for ConvectionDiffusion2D1
     * \see Stencil
     */
    void pushProlongation(const Prolongation&) {}

    /**
     * \brief does nothing for ConvectionDiffusion2D1
     * \see Stencil
     */
    void popProlongation() {}

    /**
     * \brief does nothing for ConvectionDiffusion2D1
     * \see Stencil
     */
    void pushRestriction(const Restriction&) {}

    /**
     * \brief does nothing for ConvectionDiffusion2D1
     * \see Stencil
     */
    void popRestriction() {}

    /**
     * \brief gives the max expansion of ConvectionDiffusion2D1
     * 
     * \return  1
     */
    inline size_t size() const
    {
        return 1;
    }

    /**
     * \brief returns false, because ConvectionDiffusion2D1 is not constant
     * 
     * \return  false
     */
    inline bool isConstant() const
    {
        return false;
    }
};

}

#endif /*CONVECTIONDIFFUSION2D1_H_*/
