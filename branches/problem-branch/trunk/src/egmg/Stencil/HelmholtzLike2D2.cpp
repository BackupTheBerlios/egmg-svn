/** \file HelmholtzLike2D2.cpp
 * \author Jiri Kraus
 * \brief Contains the implementation of the class HelmholtzLike2D2
 */

#include "HelmholtzLike2D2.h"

namespace mg
{

PositionArray HelmholtzLike2D2::initJx_() const
{
    const int t[]={0,-1,0,1,0};
    return PositionArray(t,5);
}

PositionArray HelmholtzLike2D2::initJy_() const
{
    const int t[]={0,0,1,0,-1};
    return PositionArray(t,5);
}

HelmholtzLike2D2::HelmholtzLike2D2(
    Precision ax,
    Precision ay,
    const DiscreteFunction& c )
    : jx_( initJx_() ), jy_( initJy_() ),
      ax_( ax ), ay_( ay ), c_( c ), indexFactor_(1)
{}

HelmholtzLike2D2::~HelmholtzLike2D2() {}

NumericArray HelmholtzLike2D2::getL(
        const Position,
        const Index sx,
        const Index sy,
        const Index nx,
        const Index ny,
        const Precision hx,
        const Precision hy,
        const Point origin) const
{
	NumericArray result( 0.0, 5 );
    result[0]=(2.0*ax_)/(hx*hx)+(2.0*ay_)/(hy*hy)+c_(indexFactor_*sx,indexFactor_*sy);
    result[1]=result[3]=(-1.0*ax_)/(hx*hx);
    result[2]=result[4]=(-1.0*ay_)/(hy*hy);
    return result;
}

PositionArray HelmholtzLike2D2::getJx(
	const Position,
	const Index,
	const Index ) const
{
    return jx_;
}

PositionArray HelmholtzLike2D2::getJy(
	const Position,
	const Index,
	const Index ) const
{
    return jy_;
}

Index HelmholtzLike2D2::size() const
{
    return 1;
}

bool HelmholtzLike2D2::isConstant() const
{
    return false;
}

void HelmholtzLike2D2::pushTransferOperators(
    const Restriction&,
    const Prolongation&,
    const Index,
    const Index)
{
    indexFactor_*=2;
}

void HelmholtzLike2D2::popTransferOperators()
{
    indexFactor_/=2;
}

}
