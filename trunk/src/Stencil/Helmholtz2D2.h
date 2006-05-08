/** \file Helmholtz2D2.h
 * \author Andr� Oeckerath
 * \brief Contains the class Helmholtz2D2
 */
#ifndef Helmholtz2D2_H_
#define Helmholtz2D2_H_


#include "Stencil.h"
#include "../Prolongation/Prolongation.h"
#include "../Restriction/Restriction.h"

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
private:
    mutable NumericArray l_;
    const PositionArray jx_;
    const PositionArray jy_;
    const Precision ax_;
    const Precision ay_;
    const Precision c_;
    //initilize jx_, makes it possible to make jx_ const
    PositionArray initJx_() const
    {
        const int t[]={0,-1,0,1,0};
        return PositionArray(t,5);
    }
    //initilize jy_, makes it possible to make jy_ const
    PositionArray initJy_() const
    {
        const int t[]={0,0,1,0,-1};
        return PositionArray(t,5);
    }
    //we don't want the autogenerated copy constructor and assignment operator
    Helmholtz2D2(const Helmholtz2D2&);
    Helmholtz2D2& operator=(const Helmholtz2D2&);
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
        Precision ax=1.0,
        Precision ay=1.0,
        Precision c=1.0) 
        : l_(5),
          jx_(initJx_()), jy_(initJy_()),
          ax_(ax), ay_(ay), c_(c)
    {}

    virtual ~Helmholtz2D2() {}

    inline Precision apply(
        const NumericArray& u,
        const Position,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny) const
    {
        return 
             (2.0*ax_*nx*nx+2.0*ay_*ny*ny+c_)*u[sy*(nx+1)+sx]
            -1.0*ax_*nx*nx*u[sy*(nx+1)+sx-1]
            -1.0*ax_*nx*nx*u[sy*(nx+1)+sx+1]
            -1.0*ay_*ny*ny*u[(sy-1)*(nx+1)+sx]
            -1.0*ay_*ny*ny*u[(sy+1)*(nx+1)+sx];
    }

    inline Precision getCenter(
        const Position,
        const Index,
        const Index,
        const Index nx,
        const Index ny) const
    {
        return 2.0*ax_*nx*nx+2.0*ay_*ny*ny+c_;
    }

    inline const NumericArray& getL(
        const Position,
        const Index,
        const Index,
        const Index nx,
        const Index ny) const
    {
        l_[0]=2.0*ax_*nx*nx+2.0*ay_*ny*ny+c_;
        l_[1]=l_[3]=-1.0*ax_*nx*nx;
        l_[2]=l_[4]=-1.0*ay_*ny*ny;
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
     * \brief does nothing for Helmholtz2D2
     * \see Stencil
     */
    virtual void pushTransferOperators(
        const Restriction&,
        const Prolongation&,
        const Index,
        const Index ) {}
    
    /**
     * \brief does nothing for Helmholtz2D2
     * \see Stencil
     */
    virtual void popTransferOperators() {}

    /**
     * \brief gives the max expansion of Helmholtz2D2
     * 
     * \return  1
     */
    inline Index size() const
    {
        return 1;
    }

    /**
     * \brief returns true, because Helmholtz2D2 is constant
     * 
     * \return  true
     */
    inline bool isConstant() const
    {
        return true;
    }
};

}

#endif /*Helmholtz2D2_H_*/
