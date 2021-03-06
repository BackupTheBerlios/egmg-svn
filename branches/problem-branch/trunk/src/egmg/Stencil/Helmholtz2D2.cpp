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

NumericArray Helmholtz2D2::getL(
        const Position,
        const Index,
        const Index,
        const Index nx,
        const Index ny,
        const Precision hx,
        const Precision hy,
        const Point origin) const
{
	NumericArray result( 0.0, 5 );
    result[0]=(2.0*ax_)/(hx*hx)+(2.0*ay_)/(hy*hy)+c_;
    result[1]=result[3]=(-1.0*ax_)/(hx*hx);
    result[2]=result[4]=(-1.0*ay_)/(hy*hy);
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
