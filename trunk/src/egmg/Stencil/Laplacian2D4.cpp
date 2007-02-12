/** \file Laplacian2D4.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class Laplacian2D4
 */

#include "Laplacian2D4.h"

namespace mg
{

std::vector<PositionArray > Laplacian2D4::initJx_()
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

std::vector<PositionArray > Laplacian2D4::initJy_()
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

Laplacian2D4::Laplacian2D4( Precision ax, Precision ay )
        : jx_( initJx_() ), jy_( initJy_() ),
          ax_( ax ), ay_( ay )
{}

Laplacian2D4::~Laplacian2D4() {}

Precision Laplacian2D4::apply(
    const NumericArray& u,
    const Position position,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny ) const
{
    switch( position )
    {
    case C:
        return 
             (30.0*ax_*nx*nx/12.0+30.0*ay_*ny*ny/12.0)*u[sy*(nx+1)+sx]
            +1.0*ax_*nx*nx/12.0*u[sy*(nx+1)+sx-2]
            -16.0*ax_*nx*nx/12.0*u[sy*(nx+1)+sx-1]
            -16.0*ax_*nx*nx/12.0*u[sy*(nx+1)+sx+1]
            +1.0*ax_*nx*nx/12.0*u[sy*(nx+1)+sx+2]
            +1.0*ay_*ny*ny/12.0*u[(sy-2)*(nx+1)+sx]
            -16.0*ay_*ny*ny/12.0*u[(sy-1)*(nx+1)+sx]
            -16.0*ay_*ny*ny/12.0*u[(sy+1)*(nx+1)+sx]
            +1.0*ay_*ny*ny/12.0*u[(sy+2)*(nx+1)+sx];
		break;
    default:
        return 
             (2.0*ax_*nx*nx+2.0*ay_*ny*ny)*u[sy*(nx+1)+sx]
            -1.0*ax_*nx*nx*u[sy*(nx+1)+sx-1]
            -1.0*ax_*nx*nx*u[sy*(nx+1)+sx+1]
            -1.0*ay_*ny*ny*u[(sy-1)*(nx+1)+sx]
            -1.0*ay_*ny*ny*u[(sy+1)*(nx+1)+sx];
    }
}

Precision Laplacian2D4::getCenter(
    const Position position,
    const Index,
    const Index,
    const Index nx,
    const Index ny ) const
{

    switch( position )
    {
    case C:
        return 30.0*ax_*nx*nx/12.0+30.0*ay_*ny*ny/12.0;
		break;
    default:
        return 2.0*ax_*nx*nx+2.0*ay_*ny*ny;
    }
}

NumericArray Laplacian2D4::getL(
    const Position position,
    const Index,
    const Index,
    const Index nx,
    const Index ny) const
{
    switch (position)
    {
    case C:
	{
		NumericArray result( 0.0, 9 );
        result[0] =30.0*ax_*nx*nx/12.0+30.0*ay_*ny*ny/12.0;
        result[1]=result[3]=-16.0*ax_*nx*nx/12.0;
        result[5]=result[7]=+1.0*ax_*nx*nx/12.0;
        result[2]=result[4]=-16.0*ay_*ny*ny/12.0;
        result[6]=result[8]=+1.0*ay_*ny*ny/12.0;
        return result;
	}
	break;
    case N:
    case E:
    case S:
    case W:
	{
		NumericArray result( 0.0, 8 );
        result[0]=2.0*ax_*nx*nx+2.0*ay_*ny*ny;
        result[1]=result[3]=-1.0*ax_*nx*nx;
        result[2]=result[4]=-1.0*ay_*ny*ny;
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
        result[0]=2.0*ax_*nx*nx+2.0*ay_*ny*ny;
        result[1]=result[3]=-1.0*ax_*nx*nx;
        result[2]=result[4]=-1.0*ay_*ny*ny;
        result[5]=result[6]=0.0;
        return result;
	}
	break;
    default:
	{
        //assert( true );
		NumericArray result( 1.0, 1);
        return result;
	}
    }
}

PositionArray Laplacian2D4::getJx(
	const Position pos,
	const Index,
	const Index ) const
{
    return jx_[pos];

}

PositionArray Laplacian2D4::getJy(
	const Position pos,
	const Index,
	const Index ) const
{
    return jy_[pos];
}


Index Laplacian2D4::size() const
{
    return 2;
}

bool Laplacian2D4::isConstant() const
{
    return true;
}

}
