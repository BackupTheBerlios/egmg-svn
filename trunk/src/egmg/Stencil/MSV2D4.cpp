/** \file MSV2D4.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class MSV2D4
 */

#include "MSV2D4.h"

namespace mg
{

PositionArray MSV2D4::initJx_() const
{
    const int t[]={0,-1,0,1,0,-1,1,-1,1};
    return PositionArray(t,9);
}

PositionArray MSV2D4::initJy_() const
{
    const int t[]={0,0,1,0,-1,-1,-1,1,1};
    return PositionArray(t,9);
}

MSV2D4::MSV2D4( Precision ax, Precision ay ) 
    : jx_( initJx_() ), jy_( initJy_() ),
	  ax_( ax ), ay_( ay )
{}

MSV2D4::~MSV2D4() {}

Precision MSV2D4::apply(
    const NumericArray& u,
    const Position,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny) const
{
    return 
         ((20.0*nx*nx+20.0*ny*ny)/12.0)*u[sy*(nx+1)+sx]
        -10.0*nx*nx/12.0*u[sy*(nx+1)+sx-1]
        -10.0*nx*nx/12.0*u[sy*(nx+1)+sx+1]
        -10.0*ny*ny/12.0*u[sy*(nx+1)+sx-(nx+1)]
        -10.0*ny*ny/12.0*u[sy*(nx+1)+sx+(nx+1)]
        +2.0*nx*nx/12.0*u[sy*(nx+1)+sx-(nx+1)]
        +2.0*nx*nx/12.0*u[sy*(nx+1)+sx+(nx+1)]
        +2.0*ny*ny/12.0*u[sy*(nx+1)+sx-1]
        +2.0*ny*ny/12.0*u[sy*(nx+1)+sx+1]
        -1.0*nx*nx/12.0*
            ( u[sy*(nx+1)+sx-(nx+1)-1]+u[sy*(nx+1)+sx-(nx+1)+1]
             +u[sy*(nx+1)+sx+(nx+1)-1]+u[sy*(nx+1)+sx+(nx+1)+1])
        -1.0*ny*ny/12.0*
            ( u[sy*(nx+1)+sx-(nx+1)-1]+u[sy*(nx+1)+sx-(nx+1)+1]
             +u[sy*(nx+1)+sx+(nx+1)-1]+u[sy*(nx+1)+sx+(nx+1)+1]);
}

Precision MSV2D4::getCenter(
    const Position,
    const Index,
    const Index,
    const Index nx,
    const Index ny) const
{
    return (20.0*nx*nx+20.0*ny*ny)/12.0;
}

NumericArray MSV2D4::getL(
    const Position,
    const Index,
    const Index,
    const Index nx,
    const Index ny) const
{
	NumericArray result( 0.0, 9 );
    result[0]=20.0*nx*nx/12+20.0*ny*ny/12;
    result[1]=result[3]=-10.0*nx*nx/12+2.0*ny*ny/12;
    result[2]=result[4]=-10.0*ny*ny/12+2.0*nx*nx/12;
    result[5]=result[6]=result[7]=result[8]=-1.0*nx*nx/12-1.0*ny*ny/12;
    return result;
}

PositionArray MSV2D4::getJx(
	const Position,
	const Index,
	const Index ) const
{
    return jx_;
}

PositionArray MSV2D4::getJy(
	const Position,
	const Index,
	const Index ) const
{
    return jy_;
}

Index MSV2D4::size() const
{
    return 1;
}

bool MSV2D4::isConstant() const
{
    return true;
}

}
