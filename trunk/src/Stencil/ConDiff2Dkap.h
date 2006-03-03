/** \file ConDiff2Dkap.h
 * \author Andr� Oeckerath
 * \brief Contains the class ConDiff2Dkap
 */
#ifndef ConDiff2Dkap_H_
#define ConDiff2Dkap_H_

#include <vector>
#include <cmath>

#include "Stencil.h"
#include "../Prolongation/Prolongation.h"
#include "../Restriction/Restriction.h"

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
private:
    mutable NumericArray lCenter_;
    mutable NumericArray lBorder_;
    mutable NumericArray lCorner_;
    const std::vector<PositionArray > jx_;
    const std::vector<PositionArray > jy_;
    const Precision epsilon_;
    const Precision factor_;
    const Precision a1_;
    const Precision a2_;
    
    std::vector<PositionArray > initJx_()
    {
        std::vector<PositionArray > jx(9);

        const int jxCenter[]={0,-1,0,1,0,-2,0,2,0};

        jx[C].resize(9);
        jx[C]=PositionArray(jxCenter,9);

        int jxBorder[]={0,-1,0,1,0,0,2,0};

        jx[W].resize(8);
        jx[W]=PositionArray(jxBorder,8);

        jxBorder[5]=-2;
        jxBorder[6]=2;
        jxBorder[7]=0;
        jx[N].resize(8);
        jx[N]=PositionArray(jxBorder,8);

        jxBorder[5]=-2;
        jxBorder[6]=0;
        jxBorder[7]=0;
        jx[E].resize(8);
        jx[E]=PositionArray(jxBorder,8);

        jxBorder[5]=-2;
        jxBorder[6]=0;
        jxBorder[7]=2;
        jx[S].resize(8);
        jx[S]=PositionArray(jxBorder,8);

        int jxCorner[]={0,-1,0,1,0,2,0};

        jx[NW].resize(7);
        jx[NW]=PositionArray(jxCorner,7);

        jxCorner[5]=-2;
        jxCorner[6]=0;
        jx[NE].resize(7);
        jx[NE]=PositionArray(jxCorner,7);

        jxCorner[5]=-2;
        jxCorner[6]=0;
        jx[SE].resize(7);
        jx[SE]=PositionArray(jxCorner,7);

        jxCorner[5]=0;
        jxCorner[6]=2;
        jx[SW].resize(7);
        jx[SW]=PositionArray(jxCorner,7);

        return jx;
    }

    std::vector<PositionArray > initJy_()
    {
        std::vector<PositionArray > jy(9);

        const int jyCenter[]={0,0,1,0,-1,0,2,0,-2};

        jy[C].resize(9);
        jy[C]=PositionArray(jyCenter,9);

        int jyBorder[]={0,0,1,0,-1,2,0,-2};

        jy[W].resize(8);
        jy[W]=PositionArray(jyBorder,8);

        jyBorder[5]=jyBorder[6]=0;
        jyBorder[7]=-2;
        jy[N].resize(8);
        jy[N]=PositionArray(jyBorder,8);

        jyBorder[5]=0;
        jyBorder[6]=2;
        jyBorder[7]=-2;
        jy[E].resize(8);
        jy[E]=PositionArray(jyBorder,8);

        jyBorder[5]=0;
        jyBorder[6]=2;
        jyBorder[7]=0;
        jy[S].resize(8);
        jy[S]=PositionArray(jyBorder,8);

        int jyCorner[]={0,0,1,0,-1,0,-2};

        jy[NW].resize(7);
        jy[NW]=PositionArray(jyCorner,7);

        jyCorner[5]=0;
        jyCorner[6]=-2;
        jy[NE].resize(7);
        jy[NE]=PositionArray(jyCorner,7);

        jyCorner[5]=0;
        jyCorner[6]=2;
        jy[SE].resize(7);
        jy[SE]=PositionArray(jyCorner,7);

        jyCorner[5]=2;
        jyCorner[6]=0;
        jy[SW].resize(7);
        jy[SW]=PositionArray(jyCorner,7);

        return jy;
    }
    //
    //we don't want the autogenerated copy constructor  and assignment operator
    ConDiff2Dkap(const ConDiff2Dkap&);
    ConDiff2Dkap& operator=(const ConDiff2Dkap&);
public:
    /**
     * \brief The constructor of a ConDiff2Dkap object
     * 
     * ConDiff2Dkap constructs a ConDiff2Dkap where \f$a_x\f$ and \f$a_y\f$
     * are given by:
     * \param[in] a_x   coefficient of the diff. operator (default 1.0)
     * \param[in] a_y   coefficient of the diff. operator (default 1.0)
     */
    explicit ConDiff2Dkap(
        Precision epsilon=1.0,
        Precision beta=0.78539816339745,
        Precision kappa=1.0/3)
            : lCenter_(9), lBorder_(8), lCorner_(7),
              jx_(initJx_()),jy_(initJy_()),
              epsilon_(epsilon), factor_((1.0-kappa)/8.0),
              a1_(cos(beta)), a2_(sin(beta))
    {}
    virtual ~ConDiff2Dkap() {}

    inline Precision apply(
        const NumericArray& u,
        const Position pos,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny) const
    {
        switch (pos)
        {
        case C:
            return 
                 (2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny
                  +factor_*(6.0*nx*fabs(a1_)+6.0*ny*fabs(a2_))
                 )*u[sy*(nx+1)+sx]
                +factor_*nx*(fabs(a1_)+a1_)*u[sy*(nx+1)+sx-2]
                +(factor_*nx*(-4.0*fabs(a1_)-2.0*a1_)-epsilon_*nx*nx
                 )*u[sy*(nx+1)+sx-1]
                +(factor_*nx*(-4.0*fabs(a1_)+2.0*a1_)-epsilon_*nx*nx
                 )*u[sy*(nx+1)+sx+1]
                +factor_*nx*(fabs(a1_)-a1_)*u[sy*(nx+1)+sx+2]
                +factor_*ny*(fabs(a2_)+a2_)*u[(sy-2)*(nx+1)+sx]
                +(factor_*ny*(-4.0*fabs(a2_)-2.0*a2_)-epsilon_*ny*ny
                 )*u[(sy-1)*(nx+1)+sx]
                +(factor_*ny*(-4.0*fabs(a2_)+2.0*a2_)-epsilon_*ny*ny
                 )*u[(sy+1)*(nx+1)+sx]
                +factor_*ny*(fabs(a2_)-a2_)*u[(sy+2)*(nx+1)+sx];
        default:
            return 
                (2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny)*u[sy*(nx+1)+sx]
                +(-1.0*epsilon_*nx*nx-a1_*nx/2)*u[sy*(nx+1)+sx-1]
                +(-1.0*epsilon_*nx*nx+a1_*nx/2)*u[sy*(nx+1)+sx+1]
                +(-1.0*epsilon_*ny*ny-a2_*ny/2)*u[(sy-1)*(nx+1)+sx]
                +(-1.0*epsilon_*ny*ny+a2_*ny/2)*u[(sy+1)*(nx+1)+sx];
        }
    }

    inline Precision getCenter(
        const Position pos,
        const Index,
        const Index,
        const Index nx,
        const Index ny) const
    {
        switch (pos)
        {
        case C:
            return 2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny
                   +factor_*(6.0*nx*fabs(a1_)+6.0*ny*fabs(a2_));
        default:
            return 2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny;
        }
    }

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
    inline const NumericArray& getL(
        const Position pos,
        const Index,
        const Index,
        const Index nx,
        const Index ny) const
    {
        switch (pos)
        {
        case C:
            lCenter_[0]=2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny
                        +factor_*(6.0*nx*fabs(a1_)+6.0*ny*fabs(a2_));
            lCenter_[1]=factor_*nx*(-4.0*fabs(a1_)-2.0*a1_)-epsilon_*nx*nx;
            lCenter_[3]=factor_*nx*(-4.0*fabs(a1_)+2.0*a1_)-epsilon_*nx*nx;
            lCenter_[2]=factor_*ny*(-4.0*fabs(a2_)+2.0*a2_)-epsilon_*ny*ny;
            lCenter_[4]=factor_*ny*(-4.0*fabs(a2_)-2.0*a2_)-epsilon_*ny*ny;
            lCenter_[5]=factor_*nx*(fabs(a1_)+a1_);
            lCenter_[7]=factor_*nx*(fabs(a1_)-a1_);
            lCenter_[6]=factor_*ny*(fabs(a2_)-a2_);
            lCenter_[8]=factor_*ny*(fabs(a2_)+a2_);
            return lCenter_;
        case W:
        case N:
        case E:
        case S:
            lBorder_[0]=2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny;
            lBorder_[1]=-1.0*epsilon_*nx*nx-a1_*nx/2;
            lBorder_[3]=-1.0*epsilon_*nx*nx+a1_*nx/2;
            lBorder_[2]=-1.0*epsilon_*ny*ny-a2_*ny/2;
            lBorder_[4]=-1.0*epsilon_*ny*ny+a2_*ny/2;
            lBorder_[5]=lBorder_[6]=lBorder_[7]=0.0;
            return lBorder_;
        case NW:
        case NE:
        case SE:
        case SW:
            lCorner_[0]=2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny;
            lCorner_[1]=-1.0*epsilon_*nx*nx-a1_*nx/2;
            lCorner_[3]=-1.0*epsilon_*nx*nx+a1_*nx/2;
            lCorner_[2]=-1.0*epsilon_*ny*ny-a2_*ny/2;
            lCorner_[4]=-1.0*epsilon_*ny*ny+a2_*ny/2;
            lCorner_[5]=lCorner_[6]=0.0;
            return lCorner_;
        default:
            lCorner_=1.0;
            return lCorner_;
        }
    }

    inline const PositionArray& getJx(const Position p) const
    {
        return jx_[p];
    }

    inline const PositionArray& getJy(const Position p) const
    {
        return jy_[p];
    }

    /**
     * \brief does nothing for ConDiff2Dkap
     * \see Stencil
     */
    void pushProlongation(const Prolongation&) {}

    /**
     * \brief does nothing for ConDiff2Dkap
     * \see Stencil
     */
    void popProlongation() {}

    /**
     * \brief does nothing for ConDiff2Dkap
     * \see Stencil
     */
    void pushRestriction(const Restriction&) {}

    /**
     * \brief does nothing for ConDiff2Dkap
     * \see Stencil
     */
    void popRestriction() {}

    /**
     * \brief gives the max expansion of ConDiff2Dkap
     * 
     * \return  2
     */
    inline Index size() const
    {
        return 2;
    }

    /**
     * \brief returns true, because ConDiff2Dkap is constant
     * 
     * \return  true
     */
    inline bool isConstant() const
    {
        return true;
    }
};

}

#endif /*ConDiff2Dkap_H_*/
