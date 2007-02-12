/** \file ConDiff2Dkap.cpp
 * \author Andre Oeckerath
 * \brief Contains the implementation of the class ConDiff2Dkap
 */

#include <cmath>

#include "ConDiff2Dkap.h"

namespace mg
{

std::vector<PositionArray > ConDiff2Dkap::initJx_()
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

std::vector<PositionArray > ConDiff2Dkap::initJy_()
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

ConDiff2Dkap::ConDiff2Dkap(
    Precision epsilon,
    Precision beta,
    Precision kappa )
        : jx_( initJx_() ), jy_( initJy_() ),
          epsilon_( epsilon ), factor_( ( 1.0 - kappa ) / 8.0 ),
          a1_( cos( beta ) ), a2_( sin( beta ) )
{}

ConDiff2Dkap::~ConDiff2Dkap() {}

Precision ConDiff2Dkap::apply(
    const NumericArray& u,
    const Position pos,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny) const
{
    switch ( pos )
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
		break;
    default:
        return 
            (2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny)*u[sy*(nx+1)+sx]
            +(-1.0*epsilon_*nx*nx-a1_*nx/2)*u[sy*(nx+1)+sx-1]
            +(-1.0*epsilon_*nx*nx+a1_*nx/2)*u[sy*(nx+1)+sx+1]
            +(-1.0*epsilon_*ny*ny-a2_*ny/2)*u[(sy-1)*(nx+1)+sx]
            +(-1.0*epsilon_*ny*ny+a2_*ny/2)*u[(sy+1)*(nx+1)+sx];
    }
}

Precision ConDiff2Dkap::getCenter(
    const Position pos,
    const Index,
    const Index,
    const Index nx,
    const Index ny) const
{
    switch ( pos )
    {
    case C:
        return 2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny
               +factor_*(6.0*nx*fabs(a1_)+6.0*ny*fabs(a2_));
		break;
    default:
        return 2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny;
    }
}

NumericArray ConDiff2Dkap::getL(
    const Position pos,
    const Index,
    const Index,
    const Index nx,
    const Index ny) const
{
    switch ( pos )
    {
    case C:
	{
		NumericArray result( 0.0, 9 );
        result[0]=2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny
                    +factor_*(6.0*nx*fabs(a1_)+6.0*ny*fabs(a2_));
        result[1]=factor_*nx*(-4.0*fabs(a1_)-2.0*a1_)-epsilon_*nx*nx;
        result[3]=factor_*nx*(-4.0*fabs(a1_)+2.0*a1_)-epsilon_*nx*nx;
        result[2]=factor_*ny*(-4.0*fabs(a2_)+2.0*a2_)-epsilon_*ny*ny;
        result[4]=factor_*ny*(-4.0*fabs(a2_)-2.0*a2_)-epsilon_*ny*ny;
        result[5]=factor_*nx*(fabs(a1_)+a1_);
        result[7]=factor_*nx*(fabs(a1_)-a1_);
        result[6]=factor_*ny*(fabs(a2_)-a2_);
        result[8]=factor_*ny*(fabs(a2_)+a2_);
        return result;
	}
	break;
    case W:
    case N:
    case E:
    case S:
	{
		NumericArray result( 0.0, 8 );
        result[0]=2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny;
        result[1]=-1.0*epsilon_*nx*nx-a1_*nx/2;
        result[3]=-1.0*epsilon_*nx*nx+a1_*nx/2;
        result[2]=-1.0*epsilon_*ny*ny-a2_*ny/2;
        result[4]=-1.0*epsilon_*ny*ny+a2_*ny/2;
        result[5]=result[6]=result[7]=0.0;
        return result;
	}
	break;
    case NW:
    case NE:
    case SE:
    case SW:
	{
		NumericArray result( 0.0, 7 );
        result[0]=2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny;
        result[1]=-1.0*epsilon_*nx*nx-a1_*nx/2;
        result[3]=-1.0*epsilon_*nx*nx+a1_*nx/2;
        result[2]=-1.0*epsilon_*ny*ny-a2_*ny/2;
        result[4]=-1.0*epsilon_*ny*ny+a2_*ny/2;
        result[5]=result[6]=0.0;
        return result;
	}
	break;
    default:
	{
        NumericArray result( 1.0, 1 );
        return result;
	}
    }
}

PositionArray ConDiff2Dkap::getJx(
	const Position p,
	const Index,
	const Index ) const
{
    return jx_[p];
}

PositionArray ConDiff2Dkap::getJy(
	const Position p,
	const Index,
	const Index ) const
{
    return jy_[p];
}

Index ConDiff2Dkap::size() const
{
    return 2;
}

bool ConDiff2Dkap::isConstant() const
{
    return true;
}

}
