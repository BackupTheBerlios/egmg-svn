/** \file Helmholtz2D2.h
 * \author Andre Oeckerath
 * \brief Contains the implementation of the class Helmholtz2D2
 */

#include "Helmholtz2D2.h"

namespace mg
{

PositionArray Helmholtz2D2::initJx_() const
{
    const int t[]={0,-1,0,1,0};
    return PositionArray(t,5);
}

PositionArray Helmholtz2D2::initJy_() const
{
    const int t[]={0,0,1,0,-1};
    return PositionArray(t,5);
}

Helmholtz2D2::Helmholtz2D2(
    Precision ax,
    Precision ay,
    Precision c ) 
    : jx_( initJx_() ), jy_( initJy_() ),
      ax_( ax ), ay_( ay ), c_( c )
{}

Helmholtz2D2::~Helmholtz2D2() {}

Precision Helmholtz2D2::apply(
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

Precision Helmholtz2D2::getCenter(
    const Position,
    const Index,
    const Index,
    const Index nx,
    const Index ny) const
{
    return 2.0*ax_*nx*nx+2.0*ay_*ny*ny+c_;
}

NumericArray Helmholtz2D2::getL(
    const Position,
    const Index,
    const Index,
    const Index nx,
    const Index ny) const
{
	NumericArray result( 0.0, 5 );
    result[0]=2.0*ax_*nx*nx+2.0*ay_*ny*ny+c_;
    result[1]=result[3]=-1.0*ax_*nx*nx;
    result[2]=result[4]=-1.0*ay_*ny*ny;
    return result;
}

PositionArray Helmholtz2D2::getJx(
	const Position,
	const Index,
	const Index ) const
{
    return jx_;
}

PositionArray Helmholtz2D2::getJy(
	const Position,
	const Index,
	const Index ) const
{
    return jy_;
}

Index Helmholtz2D2::size() const
{
    return 1;
}

bool Helmholtz2D2::isConstant() const
{
    return true;
}

}
