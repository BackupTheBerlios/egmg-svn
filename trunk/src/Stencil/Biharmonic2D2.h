/** \file Biharmonic2D2.h
 * \author Andr� Oeckerath
 * \brief Contains the class Biharmonic2D2
 */
#ifndef BIHARMONIC2D2_H_
#define BIHARMONIC2D2_H_

#include <vector>

#include "Stencil.h"
#include "../Prolongation/Prolongation.h"
#include "../Restriction/Restriction.h"

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
private:
    mutable NumericArray lCenter_;
    mutable NumericArray lBorder_;
    mutable NumericArray lCorner_;
    const std::vector<PositionArray > jx_;
    const std::vector<PositionArray > jy_;

    std::vector<PositionArray > initJx_()
    {
        std::vector<PositionArray > jx(13);

        const int jxCenter[]={0,-1,0,1,0,-2,0,2,0,-1,1,1,-1};
        jx[C].resize(13);
        jx[C]=PositionArray(jxCenter,13);

        int jxBorder[]={0,-1,0,1,0,0,2,0,-1,1,1,-1};
        jx[W].resize(12);
        jx[W]=PositionArray(jxBorder,12);

        jxBorder[5]=-2;
        jxBorder[6]=2;
        jxBorder[7]=0;
        jxBorder[8]=-1;
        jxBorder[9]=1;
        jxBorder[10]=1;
        jxBorder[11]=-1;
        jx[N].resize(12);
        jx[N]=PositionArray(jxBorder,12);

        jxBorder[5]=-2;
        jxBorder[6]=0;
        jxBorder[7]=0;
        jxBorder[8]=-1;
        jxBorder[9]=1;
        jxBorder[10]=1;
        jxBorder[11]=-1;
        jx[E].resize(12);
        jx[E]=PositionArray(jxBorder,12);

        jxBorder[5]=-2;
        jxBorder[6]=0;
        jxBorder[7]=2;
        jxBorder[8]=-1;
        jxBorder[9]=1;
        jxBorder[10]=1;
        jxBorder[11]=-1;
        jx[S].resize(12);
        jx[S]=PositionArray(jxBorder,12);

        int jxCorner[]={0,-1,0,1,0,2,0,-1,1,1,-1};
        jx[NW].resize(11);
        jx[NW]=PositionArray(jxCorner,11);

        jxCorner[5]=-2;
        jxCorner[6]=0;
        jxBorder[7]=-1;
        jxBorder[8]=1;
        jxBorder[9]=1;
        jxBorder[10]=-1;
        jx[NE].resize(11);
        jx[NE]=PositionArray(jxCorner,11);

        jxCorner[5]=-2;
        jxCorner[6]=0;
        jxBorder[7]=-1;
        jxBorder[8]=1;
        jxBorder[9]=1;
        jxBorder[10]=-1;
        jx[SE].resize(11);
        jx[SE]=PositionArray(jxCorner,11);

        jxCorner[5]=0;
        jxCorner[6]=2;
        jxBorder[7]=-1;
        jxBorder[8]=1;
        jxBorder[9]=1;
        jxBorder[10]=-1;
        jx[SW].resize(11);
        jx[SW]=PositionArray(jxCorner,11);

        return jx;
    }

    std::vector<PositionArray > initJy_()
    {
        std::vector<PositionArray > jy(13);

        const int jyCenter[]={0,0,1,0,-1,0,2,0,-2,1,1,-1,-1};

        jy[C].resize(13);
        jy[C]=PositionArray(jyCenter,13);

        int jyBorder[]={0,0,1,0,-1,2,0,-2,1,1,-1,-1};

        jy[W].resize(12);
        jy[W]=PositionArray(jyBorder,12);

        jyBorder[5]=jyBorder[6]=0;
        jyBorder[7]=-2;
        jyBorder[8]=1;
        jyBorder[9]=1;
        jyBorder[10]=-1;
        jyBorder[11]=-1;
        jy[N].resize(12);
        jy[N]=PositionArray(jyBorder,12);

        jyBorder[5]=0;
        jyBorder[6]=2;
        jyBorder[7]=-2;
        jyBorder[8]=1;
        jyBorder[9]=1;
        jyBorder[10]=-1;
        jyBorder[11]=-1;
        jy[E].resize(12);
        jy[E]=PositionArray(jyBorder,12);

        jyBorder[5]=0;
        jyBorder[6]=2;
        jyBorder[7]=0;
        jyBorder[8]=1;
        jyBorder[9]=1;
        jyBorder[10]=-1;
        jyBorder[11]=-1;
        jy[S].resize(12);
        jy[S]=PositionArray(jyBorder,12);

        int jyCorner[]={0,0,1,0,-1,0,-2,1,1,-1,-1};

        jy[NW].resize(11);
        jy[NW]=PositionArray(jyCorner,11);

        jyCorner[5]=0;
        jyCorner[6]=-2;
        jyBorder[7]=1;
        jyBorder[8]=1;
        jyBorder[9]=-1;
        jyBorder[10]=-1;
        jy[NE].resize(11);
        jy[NE]=PositionArray(jyCorner,11);

        jyCorner[5]=0;
        jyCorner[6]=2;
        jyBorder[7]=1;
        jyBorder[8]=1;
        jyBorder[9]=-1;
        jyBorder[10]=-1;
        jy[SE].resize(11);
        jy[SE]=PositionArray(jyCorner,11);

        jyCorner[5]=2;
        jyCorner[6]=0;
        jyBorder[7]=1;
        jyBorder[8]=1;
        jyBorder[9]=-1;
        jyBorder[10]=-1;
        jy[SW].resize(11);
        jy[SW]=PositionArray(jyCorner,11);

        return jy;
    }
    //we don't want the autogenerated copy constructor and assignment operator
    Biharmonic2D2(const Biharmonic2D2&);
    Biharmonic2D2& operator=(const Biharmonic2D2&);
public:
    /**
     * \brief The constructor of a Biharmonic2D2 object
     * 
     * Biharmonic2D2 constructs a Biharmonic2D2 where \f$a_x\f$ and \f$a_y\f$
     * are given by:
     * \param[in] a_x   coefficient of the diff. operator (default 1.0)
     * \param[in] a_y   coefficient of the diff. operator (default 1.0)
     */
    explicit Biharmonic2D2()
            : lCenter_(13), lBorder_(12), lCorner_(11),
              jx_(initJx_()), jy_(initJy_()) {}
    virtual ~Biharmonic2D2() {}

    inline Precision apply(
        const NumericArray& u,
        const Position pos,
        const size_t sx,
        const size_t sy,
        const size_t nx,
        const size_t ny) const
    {
        switch(pos)
        {
        case C: return
            (6.0*nx*nx*nx*nx+6.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx]
            +1.0*nx*nx*nx*nx*u[sy*(nx+1)+sx-2]
            -(4.0*nx*nx*nx*nx+4.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx-1]
            -(4.0*nx*nx*nx*nx+4.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx+1]
            +1.0*nx*nx*nx*nx*u[sy*(nx+1)+sx+2]
            +1.0*ny*ny*ny*ny*u[(sy-2)*(nx+1)+sx]
            -(4.0*ny*ny*ny*ny+4.0*nx*nx*ny*ny)*u[(sy-1)*(nx+1)+sx]
            -(4.0*ny*ny*ny*ny+4.0*nx*nx*ny*ny)*u[(sy+1)*(nx+1)+sx]
            +1.0*ny*ny*ny*ny*u[(sy+2)*(nx+1)+sx]
            +2.0*nx*nx*ny*ny*u[(sy+1)*(nx+1)+sx-1]
            +2.0*nx*nx*ny*ny*u[(sy+1)*(nx+1)+sx+1]
            +2.0*nx*nx*ny*ny*u[(sy-1)*(nx+1)+sx+1]
            +2.0*nx*nx*ny*ny*u[(sy-1)*(nx+1)+sx-1];
        case W: return
            (7.0*nx*nx*nx*nx+6.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx]
            -(4.0*nx*nx*nx*nx+4.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx-1]
            -(4.0*nx*nx*nx*nx+4.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx+1]
            +1.0*nx*nx*nx*nx*u[sy*(nx+1)+sx+2]
            +1.0*ny*ny*ny*ny*u[(sy-2)*(nx+1)+sx]
            -(4.0*ny*ny*ny*ny+4.0*nx*nx*ny*ny)*u[(sy-1)*(nx+1)+sx]
            -(4.0*ny*ny*ny*ny+4.0*nx*nx*ny*ny)*u[(sy+1)*(nx+1)+sx]
            +1.0*ny*ny*ny*ny*u[(sy+2)*(nx+1)+sx]
            +2.0*nx*nx*ny*ny*u[(sy+1)*(nx+1)+sx-1]
            +2.0*nx*nx*ny*ny*u[(sy+1)*(nx+1)+sx+1]
            +2.0*nx*nx*ny*ny*u[(sy-1)*(nx+1)+sx+1]
            +2.0*nx*nx*ny*ny*u[(sy-1)*(nx+1)+sx-1];
        case NW: return
             (7.0*nx*nx*nx*nx+7.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx]
            -(4.0*nx*nx*nx*nx+4.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx-1]
            -(4.0*nx*nx*nx*nx+4.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx+1]
            +1.0*nx*nx*nx*nx*u[sy*(nx+1)+sx+2]
            +1.0*ny*ny*ny*ny*u[(sy-2)*(nx+1)+sx]
            -(4.0*ny*ny*ny*ny+4.0*nx*nx*ny*ny)*u[(sy-1)*(nx+1)+sx]
            -(4.0*ny*ny*ny*ny+4.0*nx*nx*ny*ny)*u[(sy+1)*(nx+1)+sx]
            +2.0*nx*nx*ny*ny*u[(sy+1)*(nx+1)+sx-1]
            +2.0*nx*nx*ny*ny*u[(sy+1)*(nx+1)+sx+1]
            +2.0*nx*nx*ny*ny*u[(sy-1)*(nx+1)+sx+1]
            +2.0*nx*nx*ny*ny*u[(sy-1)*(nx+1)+sx-1];
        case N: return
            (6.0*nx*nx*nx*nx+7.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx]
            +1.0*nx*nx*nx*nx*u[sy*(nx+1)+sx-2]
            -(4.0*nx*nx*nx*nx+4.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx-1]
            -(4.0*nx*nx*nx*nx+4.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx+1]
            +1.0*nx*nx*nx*nx*u[sy*(nx+1)+sx+2]
            +1.0*ny*ny*ny*ny*u[(sy-2)*(nx+1)+sx]
            -(4.0*ny*ny*ny*ny+4.0*nx*nx*ny*ny)*u[(sy-1)*(nx+1)+sx]
            -(4.0*ny*ny*ny*ny+4.0*nx*nx*ny*ny)*u[(sy+1)*(nx+1)+sx]
            +2.0*nx*nx*ny*ny*u[(sy+1)*(nx+1)+sx-1]
            +2.0*nx*nx*ny*ny*u[(sy+1)*(nx+1)+sx+1]
            +2.0*nx*nx*ny*ny*u[(sy-1)*(nx+1)+sx+1]
            +2.0*nx*nx*ny*ny*u[(sy-1)*(nx+1)+sx-1];
        case E: return
            (7.0*nx*nx*nx*nx+6.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx]
            +1.0*nx*nx*nx*nx*u[sy*(nx+1)+sx-2]
            -(4.0*nx*nx*nx*nx+4.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx-1]
            -(4.0*nx*nx*nx*nx+4.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx+1]
            +1.0*ny*ny*ny*ny*u[(sy-2)*(nx+1)+sx]
            -(4.0*ny*ny*ny*ny+4.0*nx*nx*ny*ny)*u[(sy-1)*(nx+1)+sx]
            -(4.0*ny*ny*ny*ny+4.0*nx*nx*ny*ny)*u[(sy+1)*(nx+1)+sx]
            +1.0*ny*ny*ny*ny*u[(sy+2)*(nx+1)+sx]
            +2.0*nx*nx*ny*ny*u[(sy+1)*(nx+1)+sx-1]
            +2.0*nx*nx*ny*ny*u[(sy+1)*(nx+1)+sx+1]
            +2.0*nx*nx*ny*ny*u[(sy-1)*(nx+1)+sx+1]
            +2.0*nx*nx*ny*ny*u[(sy-1)*(nx+1)+sx-1];
        case NE: return
            (7.0*nx*nx*nx*nx+7.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx]
            +1.0*nx*nx*nx*nx*u[sy*(nx+1)+sx-2]
            -(4.0*nx*nx*nx*nx+4.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx-1]
            -(4.0*nx*nx*nx*nx+4.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx+1]
            +1.0*ny*ny*ny*ny*u[(sy-2)*(nx+1)+sx]
            -(4.0*ny*ny*ny*ny+4.0*nx*nx*ny*ny)*u[(sy-1)*(nx+1)+sx]
            -(4.0*ny*ny*ny*ny+4.0*nx*nx*ny*ny)*u[(sy+1)*(nx+1)+sx]
            +2.0*nx*nx*ny*ny*u[(sy+1)*(nx+1)+sx-1]
            +2.0*nx*nx*ny*ny*u[(sy+1)*(nx+1)+sx+1]
            +2.0*nx*nx*ny*ny*u[(sy-1)*(nx+1)+sx+1]
            +2.0*nx*nx*ny*ny*u[(sy-1)*(nx+1)+sx-1];
        case SE: return
            (7.0*nx*nx*nx*nx+7.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx]
            +1.0*nx*nx*nx*nx*u[sy*(nx+1)+sx-2]
            -(4.0*nx*nx*nx*nx+4.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx-1]
            -(4.0*nx*nx*nx*nx+4.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx+1]
            -(4.0*ny*ny*ny*ny+4.0*nx*nx*ny*ny)*u[(sy-1)*(nx+1)+sx]
            -(4.0*ny*ny*ny*ny+4.0*nx*nx*ny*ny)*u[(sy+1)*(nx+1)+sx]
            +1.0*ny*ny*ny*ny*u[(sy+2)*(nx+1)+sx]
            +2.0*nx*nx*ny*ny*u[(sy+1)*(nx+1)+sx-1]
            +2.0*nx*nx*ny*ny*u[(sy+1)*(nx+1)+sx+1]
            +2.0*nx*nx*ny*ny*u[(sy-1)*(nx+1)+sx+1]
            +2.0*nx*nx*ny*ny*u[(sy-1)*(nx+1)+sx-1];
        case S: return
            (6.0*nx*nx*nx*nx+7.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx]
            +1.0*nx*nx*nx*nx*u[sy*(nx+1)+sx-2]
            -(4.0*nx*nx*nx*nx+4.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx-1]
            -(4.0*nx*nx*nx*nx+4.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx+1]
            +1.0*nx*nx*nx*nx*u[sy*(nx+1)+sx+2]
            -(4.0*ny*ny*ny*ny+4.0*nx*nx*ny*ny)*u[(sy-1)*(nx+1)+sx]
            -(4.0*ny*ny*ny*ny+4.0*nx*nx*ny*ny)*u[(sy+1)*(nx+1)+sx]
            +1.0*ny*ny*ny*ny*u[(sy+2)*(nx+1)+sx]
            +2.0*nx*nx*ny*ny*u[(sy+1)*(nx+1)+sx-1]
            +2.0*nx*nx*ny*ny*u[(sy+1)*(nx+1)+sx+1]
            +2.0*nx*nx*ny*ny*u[(sy-1)*(nx+1)+sx+1]
            +2.0*nx*nx*ny*ny*u[(sy-1)*(nx+1)+sx-1];
        case SW: return
            (7.0*nx*nx*nx*nx+7.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx]
            -(4.0*nx*nx*nx*nx+4.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx-1]
            -(4.0*nx*nx*nx*nx+4.0*nx*nx*ny*ny)*u[sy*(nx+1)+sx+1]
            +1.0*nx*nx*nx*nx*u[sy*(nx+1)+sx+2]
            -(4.0*ny*ny*ny*ny+4.0*nx*nx*ny*ny)*u[(sy-1)*(nx+1)+sx]
            -(4.0*ny*ny*ny*ny+4.0*nx*nx*ny*ny)*u[(sy+1)*(nx+1)+sx]
            +1.0*ny*ny*ny*ny*u[(sy+2)*(nx+1)+sx]
            +2.0*nx*nx*ny*ny*u[(sy+1)*(nx+1)+sx-1]
            +2.0*nx*nx*ny*ny*u[(sy+1)*(nx+1)+sx+1]
            +2.0*nx*nx*ny*ny*u[(sy-1)*(nx+1)+sx+1]
            +2.0*nx*nx*ny*ny*u[(sy-1)*(nx+1)+sx-1];
        default: return 1.0;
        }
    }

    inline Precision getCenter(
        const Position pos,
        const size_t,
        const size_t,
        const size_t nx,
        const size_t ny) const
    {
        switch (pos)
        {
        case C:
            return (6.0*nx*nx*nx*nx+6.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny);
        case W:
        case E:
            return (7.0*nx*nx*nx*nx+6.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny);
        case N:
        case S:
            return (6.0*nx*nx*nx*nx+7.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny);   
        case NW:
        case SW:
        case NE:
        case SE:
            return (7.0*nx*nx*nx*nx+7.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny);  
        default:
            return 1.0; 
        }
    }

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
     * \param[in]       the x coordinate of the center element (not used)
     * \param[in]       the y coordinate of the center element (not used)
     * \param[in] nx    the step size in x direction
     * \param[in] ny    the step size in y direction
     * \return          the coefficients of Biharmonic2D4
     */
    inline const NumericArray& getL(
        const Position pos,
        const size_t,
        const size_t,
        const size_t nx,
        const size_t ny) const
    {
        switch (pos)
        {
        case C:
            lCenter_[0]=6.0*nx*nx*nx*nx+6.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny;
            lCenter_[1]=lCenter_[3]=-4.0*nx*nx*nx*nx-4.0*nx*nx*ny*ny;
            lCenter_[5]=lCenter_[7]=+1.0*nx*nx*nx*nx;
            lCenter_[2]=lCenter_[4]=-4.0*ny*ny*ny*ny-4.0*nx*nx*ny*ny;
            lCenter_[6]=lCenter_[8]=+1.0*ny*ny*ny*ny;
            lCenter_[9]=lCenter_[10]=lCenter_[11]=lCenter_[12]=+2.0*nx*nx*ny*ny;
            return lCenter_;
        case W:
            lBorder_[0]=7.0*nx*nx*nx*nx+6.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny;
            lBorder_[1]=lBorder_[3]=-4.0*nx*nx*nx*nx-4.0*nx*nx*ny*ny;
            lBorder_[2]=lBorder_[4]=-4.0*ny*ny*ny*ny-4.0*nx*nx*ny*ny;
            lBorder_[5]=lBorder_[7]=+1.0*ny*ny*ny*ny;
            lBorder_[6]=+1.0*nx*nx*nx*nx;
            lBorder_[8]=lBorder_[9]=+2.0*nx*nx*ny*ny;
            lBorder_[10]=lBorder_[11]=+2.0*nx*nx*ny*ny;
            return lBorder_;
        case NW:
            lCorner_[0]=7.0*nx*nx*nx*nx+7.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny;
            lCorner_[1]=lCorner_[3]=-4.0*nx*nx*nx*nx-4.0*nx*nx*ny*ny;
            lCorner_[2]=lCorner_[4]=-4.0*ny*ny*ny*ny-4.0*nx*nx*ny*ny;
            lCorner_[5]=+1.0*nx*nx*nx*nx;
            lCorner_[6]=+1.0*ny*ny*ny*ny;
            lCorner_[7]=lCorner_[8]=+2.0*nx*nx*ny*ny;
            lCorner_[9]=lCorner_[10]=+2.0*nx*nx*ny*ny;
            return lCorner_;
        case N:
            lBorder_[0]=6.0*nx*nx*nx*nx+7.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny;
            lBorder_[1]=lBorder_[3]=-4.0*nx*nx*nx*nx-4.0*nx*nx*ny*ny;
            lBorder_[2]=lBorder_[4]=-4.0*ny*ny*ny*ny-4.0*nx*nx*ny*ny;
            lBorder_[5]=lBorder_[6]=+1.0*nx*nx*nx*nx;
            lBorder_[7]=+1.0*ny*ny*ny*ny;
            lBorder_[8]=lBorder_[9]=+2.0*nx*nx*ny*ny;
            lBorder_[10]=lBorder_[11]=+2.0*nx*nx*ny*ny;
            return lBorder_;
        case NE:
            lCorner_[0]=7.0*nx*nx*nx*nx+7.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny;
            lCorner_[1]=lCorner_[3]=-4.0*nx*nx*nx*nx-4.0*nx*nx*ny*ny;
            lCorner_[2]=lCorner_[4]=-4.0*ny*ny*ny*ny-4.0*nx*nx*ny*ny;
            lCorner_[5]=+1.0*nx*nx*nx*nx;
            lCorner_[6]=+1.0*ny*ny*ny*ny;
            lCorner_[7]=lCorner_[8]=+2.0*nx*nx*ny*ny;
            lCorner_[9]=lCorner_[10]=+2.0*nx*nx*ny*ny;
            return lCorner_;
        case E:
            lBorder_[0]=7.0*nx*nx*nx*nx+6.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny;
            lBorder_[1]=lBorder_[3]=-4.0*nx*nx*nx*nx-4.0*nx*nx*ny*ny;
            lBorder_[2]=lBorder_[4]=-4.0*ny*ny*ny*ny-4.0*nx*nx*ny*ny;
            lBorder_[6]=lBorder_[7]=+1.0*ny*ny*ny*ny;
            lBorder_[5]=+1.0*nx*nx*nx*nx;
            lBorder_[8]=lBorder_[9]=+2.0*nx*nx*ny*ny;
            lBorder_[10]=lBorder_[11]=+2.0*nx*nx*ny*ny;
            return lBorder_;
        case SE:
            lCorner_[0]=7.0*nx*nx*nx*nx+7.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny;
            lCorner_[1]=lCorner_[3]=-4.0*nx*nx*nx*nx-4.0*nx*nx*ny*ny;
            lCorner_[2]=lCorner_[4]=-4.0*ny*ny*ny*ny-4.0*nx*nx*ny*ny;
            lCorner_[5]=+1.0*nx*nx*nx*nx;
            lCorner_[6]=+1.0*ny*ny*ny*ny;
            lCorner_[7]=lCorner_[8]=+2.0*nx*nx*ny*ny;
            lCorner_[9]=lCorner_[10]=+2.0*nx*nx*ny*ny;
            return lCorner_;
        case S:
            lBorder_[0]=7.0*nx*nx*nx*nx+6.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny;
            lBorder_[1]=lBorder_[3]=-4.0*nx*nx*nx*nx-4.0*nx*nx*ny*ny;
            lBorder_[2]=lBorder_[4]=-4.0*ny*ny*ny*ny-4.0*nx*nx*ny*ny;
            lBorder_[5]=lBorder_[7]=+1.0*nx*nx*nx*nx;
            lBorder_[6]=+1.0*ny*ny*ny*ny;
            lBorder_[8]=lBorder_[9]=+2.0*nx*nx*ny*ny;
            lBorder_[10]=lBorder_[11]=+2.0*nx*nx*ny*ny;
            return lBorder_;
        case SW:
            lCorner_[0]=7.0*nx*nx*nx*nx+7.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny;
            lCorner_[1]=lCorner_[3]=-4.0*nx*nx*nx*nx-4.0*nx*nx*ny*ny;
            lCorner_[2]=lCorner_[4]=-4.0*ny*ny*ny*ny-4.0*nx*nx*ny*ny;
            lCorner_[6]=+1.0*nx*nx*nx*nx;
            lCorner_[5]=+1.0*ny*ny*ny*ny;
            lCorner_[7]=lCorner_[8]=+2.0*nx*nx*ny*ny;
            lCorner_[9]=lCorner_[10]=+2.0*nx*nx*ny*ny;
            return lCorner_;
        default:
            lCorner_=1.0;
            return lCorner_;
        }
    }

    inline const PositionArray& getJx(const Position pos) const
    {
        return jx_[pos];
    }

    inline const PositionArray& getJy(const Position pos) const
    {
        return jy_[pos];
    }

    /**
     * \brief does nothing for Biharmonic2D2
     * \see Stencil
     */
    void pushProlongation(const Prolongation&) {}

    /**
     * \brief does nothing for Biharmonic2D2
     * \see Stencil
     */
    void popProlongation() {}

    /**
     * \brief does nothing for Biharmonic2D2
     * \see Stencil
     */
    void pushRestriction(const Restriction&) {}

    /**
     * \brief does nothing for Biharmonic2D2
     * \see Stencil
     */
    void popRestriction() {}

    /**
     * \brief gives the max expansion of Biharmonic2D2
     * 
     * \return  2
     */
    inline size_t size() const
    {
        return 2;
    }

    /**
     * \brief returns true, because Biharmonic2D2 is constant
     * 
     * \return  true
     */
    inline bool isConstant() const
    {
        return true;
    }
};

}

#endif /*BIHARMONIC2D2_H_*/
