/** \file Laplacian2D2.h
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the class Laplacian2D2
 */
#ifndef LAPLACIAN2D2_H_
#define LAPLACIAN2D2_H_

#include "Stencil.h"
#include "../Prolongation/Prolongation.h"
#include "../Restriction/Restriction.h"


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
private:
    mutable NumericArray l_;
    const PositionArray jx_;
    const PositionArray jy_;
    const Precision ax_;
    const Precision ay_;
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
    Laplacian2D2( const Laplacian2D2& rhs );
    Laplacian2D2& operator=(const Laplacian2D2&);
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
        Precision ax=1.0,
        Precision ay=1.0) 
        : l_(5),
          jx_(initJx_()), jy_(initJy_()),
          ax_(ax), ay_(ay) {}

    virtual ~Laplacian2D2() {}

    inline Precision apply(
        const NumericArray& u,
        const Position,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny) const
    {
        return 
             (2.0*ax_*nx*nx+2.0*ay_*ny*ny)*u[sy*(nx+1)+sx]
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
        return 2.0*ax_*nx*nx+2.0*ay_*ny*ny;
    }

    inline const NumericArray& getL(
        const Position,
        const Index,
        const Index,
        const Index nx,
        const Index ny) const
    {
        l_[0]=2.0*ax_*nx*nx+2.0*ay_*ny*ny;
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
     * \brief does nothing for Laplacian2D2
     * \see Stencil
     */
    virtual void pushTransferOperators(
        const Restriction&,
        const Prolongation&,
        const Index,
        const Index ) {}
    
    /**
     * \brief does nothing for Laplacian2D2
     * \see Stencil
     */
    virtual void popTransferOperators() {}

    /**
     * \brief gives the max expansion of Laplacian2D2
     * 
     * \return  1
     */
    inline Index size() const
    {
        return 1;
    }

    /**
     * \brief returns true, because Laplacian2D2 is constant
     * 
     * \return  true
     */
    inline bool isConstant() const
    {
        return true;
    }
};

}

#endif /*LAPLACIAN2D2_H_*/
