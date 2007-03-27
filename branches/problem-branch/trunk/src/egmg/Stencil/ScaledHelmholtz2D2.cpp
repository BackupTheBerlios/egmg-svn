/** \file ScaledHelmholtz2D2.cpp
 * \author Jiri Kraus
 * \brief Contains the implementation of the class Helmholtz2D2
 */

#include "ScaledHelmholtz2D2.h"

namespace mg
{

PositionArray ScaledHelmholtz2D2::initJx_() const
{
    const int t[]={0,-1,0,1,0};
    return PositionArray(t,5);
}

PositionArray ScaledHelmholtz2D2::initJy_() const
{
    const int t[]={0,0,1,0,-1};
    return PositionArray(t,5);
}

NumericArray ScaledHelmholtz2D2::initL_(
    const Precision ax,
    const Precision ay,
    const Precision c) const
{
    const Precision t[]={2.0*ax+2.0*ay+c,ax,ay,ax,ay};
    return NumericArray(t,5);
}

ScaledHelmholtz2D2::ScaledHelmholtz2D2(
    Precision ax,
    Precision ay,
    Precision c )
    : jx_( initJx_() ), jy_( initJy_() ), l_( initL_(ax,ay,c) ),
      ax_( ax ), ay_( ay ), c_( c )
{}

ScaledHelmholtz2D2::~ScaledHelmholtz2D2() {}

NumericArray ScaledHelmholtz2D2::getL(
        const Position,
        const Index,
        const Index,
        const Index,
        const Index,
        const Precision,
        const Precision,
        const Point) const
{
    return l_;
}

PositionArray ScaledHelmholtz2D2::getJx(
	const Position,
	const Index,
	const Index ) const
{
    return jx_;
}

PositionArray ScaledHelmholtz2D2::getJy(
	const Position,
	const Index,
	const Index ) const
{
    return jy_;
}

Index ScaledHelmholtz2D2::size() const
{
    return 1;
}

bool ScaledHelmholtz2D2::isConstant() const
{
    return true;
}

}
