/** \file Laplacian2D2.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class Laplacian2D2
 */

#include "Laplacian2D2.h"

namespace mg
{

PositionArray Laplacian2D2::initJx_() const
{
    const int t[]={0,-1,0,1,0};
    return PositionArray(t,5);
}
PositionArray Laplacian2D2::initJy_() const
{
    const int t[]={0,0,1,0,-1};
    return PositionArray(t,5);
}

Laplacian2D2::Laplacian2D2(
    Precision ax,
    Precision ay ) 
    : jx_( initJx_() ), jy_( initJy_() ),
      ax_( ax ), ay_( ay ) {}

Laplacian2D2::~Laplacian2D2() {}

Precision Laplacian2D2::apply(
    const DiscreteFunction& u,
    const Position,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny) const
{
    return 
         (2.0*ax_*nx*nx+2.0*ay_*ny*ny)*u(sx,sy)
        -1.0*ax_*nx*nx*u(sx-1,sy)
        -1.0*ax_*nx*nx*u(sx+1,sy)
        -1.0*ay_*ny*ny*u(sx,sy-1)
        -1.0*ay_*ny*ny*u(sx,sy+1);
}

Precision Laplacian2D2::getCenter(
    const Position,
    const Index,
    const Index,
    const Index nx,
    const Index ny) const
{
    return 2.0*ax_*nx*nx+2.0*ay_*ny*ny;
}

NumericArray Laplacian2D2::getL(
    const Position,
    const Index,
    const Index,
    const Index nx,
    const Index ny) const
{
	NumericArray result( 0.0, 5 );
    result[0]=2.0*ax_*nx*nx+2.0*ay_*ny*ny;
    result[1]=result[3]=-1.0*ax_*nx*nx;
    result[2]=result[4]=-1.0*ay_*ny*ny;
    return result;
}

PositionArray Laplacian2D2::getJx(
	const Position,
	const Index,
	const Index ) const
{
    return jx_;
}

PositionArray Laplacian2D2::getJy(
	const Position,
	const Index,
	const Index ) const
{
    return jy_;
}

Index Laplacian2D2::size() const
{
    return 1;
}

bool Laplacian2D2::isConstant() const
{
    return true;
}

}
