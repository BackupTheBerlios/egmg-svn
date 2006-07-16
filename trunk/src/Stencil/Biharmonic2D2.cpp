/** \file Biharmonic2D2.cpp
 * \author André Oeckerath
 * \brief Contains implementation of the class Biharmonic2D2
 */
#include "Biharmonic2D2.h"

namespace mg
{

std::vector<PositionArray > Biharmonic2D2::initJx_()
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

std::vector<PositionArray > Biharmonic2D2::initJy_()
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

Biharmonic2D2::Biharmonic2D2()
	: jx_( initJx_() ), jy_( initJy_() ) {}

Biharmonic2D2::~Biharmonic2D2() {}

Precision Biharmonic2D2::apply(
    const NumericArray& u,
    const Position pos,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny ) const
{
    switch( pos )
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

Precision Biharmonic2D2::getCenter(
    const Position pos,
    const Index,
    const Index,
    const Index nx,
    const Index ny ) const
{
    switch ( pos )
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

NumericArray Biharmonic2D2::getL(
    const Position pos,
    const Index,
    const Index,
    const Index nx,
    const Index ny ) const
{
    switch ( pos )
    {
    case C:
	{
		NumericArray result( 0.0, 13 );
        result[0]=6.0*nx*nx*nx*nx+6.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny;
        result[1]=result[3]=-4.0*nx*nx*nx*nx-4.0*nx*nx*ny*ny;
        result[5]=result[7]=+1.0*nx*nx*nx*nx;
        result[2]=result[4]=-4.0*ny*ny*ny*ny-4.0*nx*nx*ny*ny;
        result[6]=result[8]=+1.0*ny*ny*ny*ny;
        result[9]=result[10]=result[11]=result[12]=+2.0*nx*nx*ny*ny;
        return result;
	}
	break;
    case W:
	{
		NumericArray result( 0.0, 12 );
        result[0]=7.0*nx*nx*nx*nx+6.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny;
        result[1]=result[3]=-4.0*nx*nx*nx*nx-4.0*nx*nx*ny*ny;
        result[2]=result[4]=-4.0*ny*ny*ny*ny-4.0*nx*nx*ny*ny;
        result[5]=result[7]=+1.0*ny*ny*ny*ny;
        result[6]=+1.0*nx*nx*nx*nx;
        result[8]=result[9]=+2.0*nx*nx*ny*ny;
        result[10]=result[11]=+2.0*nx*nx*ny*ny;
        return result;
	}
	break;
    case NW:
	{
		NumericArray result( 0.0, 11 );
        result[0]=7.0*nx*nx*nx*nx+7.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny;
        result[1]=result[3]=-4.0*nx*nx*nx*nx-4.0*nx*nx*ny*ny;
        result[2]=result[4]=-4.0*ny*ny*ny*ny-4.0*nx*nx*ny*ny;
        result[5]=+1.0*nx*nx*nx*nx;
        result[6]=+1.0*ny*ny*ny*ny;
        result[7]=result[8]=+2.0*nx*nx*ny*ny;
        result[9]=result[10]=+2.0*nx*nx*ny*ny;
        return result;
	}
	break;
    case N:
	{
		NumericArray result( 0.0, 12 );
        result[0]=6.0*nx*nx*nx*nx+7.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny;
        result[1]=result[3]=-4.0*nx*nx*nx*nx-4.0*nx*nx*ny*ny;
        result[2]=result[4]=-4.0*ny*ny*ny*ny-4.0*nx*nx*ny*ny;
        result[5]=result[6]=+1.0*nx*nx*nx*nx;
        result[7]=+1.0*ny*ny*ny*ny;
        result[8]=result[9]=+2.0*nx*nx*ny*ny;
        result[10]=result[11]=+2.0*nx*nx*ny*ny;
        return result;
	}
	break;
    case NE:
	{
		NumericArray result( 0.0, 11 );
        result[0]=7.0*nx*nx*nx*nx+7.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny;
        result[1]=result[3]=-4.0*nx*nx*nx*nx-4.0*nx*nx*ny*ny;
        result[2]=result[4]=-4.0*ny*ny*ny*ny-4.0*nx*nx*ny*ny;
        result[5]=+1.0*nx*nx*nx*nx;
        result[6]=+1.0*ny*ny*ny*ny;
        result[7]=result[8]=+2.0*nx*nx*ny*ny;
        result[9]=result[10]=+2.0*nx*nx*ny*ny;
        return result;
	}
	break;
    case E:
	{
		NumericArray result( 0.0, 12 );
        result[0]=7.0*nx*nx*nx*nx+6.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny;
        result[1]=result[3]=-4.0*nx*nx*nx*nx-4.0*nx*nx*ny*ny;
        result[2]=result[4]=-4.0*ny*ny*ny*ny-4.0*nx*nx*ny*ny;
        result[6]=result[7]=+1.0*ny*ny*ny*ny;
        result[5]=+1.0*nx*nx*nx*nx;
        result[8]=result[9]=+2.0*nx*nx*ny*ny;
        result[10]=result[11]=+2.0*nx*nx*ny*ny;
        return result;
	}
	break;
    case SE:
	{
		NumericArray result( 0.0, 11 );
        result[0]=7.0*nx*nx*nx*nx+7.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny;
        result[1]=result[3]=-4.0*nx*nx*nx*nx-4.0*nx*nx*ny*ny;
        result[2]=result[4]=-4.0*ny*ny*ny*ny-4.0*nx*nx*ny*ny;
        result[5]=+1.0*nx*nx*nx*nx;
        result[6]=+1.0*ny*ny*ny*ny;
        result[7]=result[8]=+2.0*nx*nx*ny*ny;
        result[9]=result[10]=+2.0*nx*nx*ny*ny;
        return result;
	}
	break;
    case S:
	{
		NumericArray result( 0.0, 12 );
        result[0]=7.0*nx*nx*nx*nx+6.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny;
        result[1]=result[3]=-4.0*nx*nx*nx*nx-4.0*nx*nx*ny*ny;
        result[2]=result[4]=-4.0*ny*ny*ny*ny-4.0*nx*nx*ny*ny;
        result[5]=result[7]=+1.0*nx*nx*nx*nx;
        result[6]=+1.0*ny*ny*ny*ny;
        result[8]=result[9]=+2.0*nx*nx*ny*ny;
        result[10]=result[11]=+2.0*nx*nx*ny*ny;
        return result;
	}
	break;
    case SW:
	{
		NumericArray result( 0.0, 11 );
        result[0]=7.0*nx*nx*nx*nx+7.0*ny*ny*ny*ny+8.0*nx*nx*ny*ny;
        result[1]=result[3]=-4.0*nx*nx*nx*nx-4.0*nx*nx*ny*ny;
        result[2]=result[4]=-4.0*ny*ny*ny*ny-4.0*nx*nx*ny*ny;
        result[6]=+1.0*nx*nx*nx*nx;
        result[5]=+1.0*ny*ny*ny*ny;
        result[7]=result[8]=+2.0*nx*nx*ny*ny;
        result[9]=result[10]=+2.0*nx*nx*ny*ny;
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

PositionArray Biharmonic2D2::getJx(
	const Position pos,
	const Index,
	const Index ) const
{
    return jx_[pos];
}

PositionArray Biharmonic2D2::getJy(
	const Position pos,
	const Index,
	const Index ) const
{
    return jy_[pos];
}

Index Biharmonic2D2::size() const
{
    return 2;
}

bool Biharmonic2D2::isConstant() const
{
    return true;
}

}
