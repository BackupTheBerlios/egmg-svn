/** \file ConDiff2D1.cpp
 * \author Andre Oeckerath
 * \brief Contains the implementation of the class ConDiff2D1
 */

#include <cmath>

#include "ConDiff2D1.h"

namespace mg
{

PositionArray ConDiff2D1::initJx_() const
{
    const int t[] = {0,-1,0,1,0};
    return PositionArray(t,5);
}

PositionArray ConDiff2D1::initJy_() const
{
    const int t[] = {0,0,1,0,-1};
    return PositionArray(t,5);
}

ConDiff2D1::ConDiff2D1(
    Precision epsilon,
    Precision beta ) 
        : jx_( initJx_() ), jy_( initJy_() ),
          epsilon_( epsilon ),
          a1_( std::cos( beta ) ), a2_( std::sin( beta ) )
{}

ConDiff2D1::~ConDiff2D1() {}

Precision ConDiff2D1::apply(
    const NumericArray& u,
    const Position,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny) const
{
    return 
        (2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny
         +fabs(a1_)*nx+fabs(a2_)*ny
        )*u[sy*(nx+1)+sx]
        +(-1.0*epsilon_*nx*nx+(-a1_-fabs(a1_))*nx/2)*u[sy*(nx+1)+sx-1]
        +(-1.0*epsilon_*nx*nx+(a1_-fabs(a1_))*nx/2)*u[sy*(nx+1)+sx+1]
        +(-1.0*epsilon_*ny*ny+(-a2_-fabs(a2_))*ny/2)*u[(sy-1)*(nx+1)+sx]
        +(-1.0*epsilon_*ny*ny+(a2_-fabs(a2_))*ny/2)*u[(sy+1)*(nx+1)+sx];
}

Precision ConDiff2D1::getCenter(
    const Position,
    const Index,
    const Index,
    const Index nx,
    const Index ny) const
{
    return 
        (2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny+fabs(a1_)*nx+fabs(a2_)*ny);
}

NumericArray ConDiff2D1::getL(
    const Position,
    const Index,
    const Index,
    const Index nx,
    const Index ny) const
{
	NumericArray result( 0.0, 5 );
    result[0] = 2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny+fabs(a1_)*nx+fabs(a2_)*ny;
    result[1] = -1.0*epsilon_*nx*nx+(-a1_-fabs(a1_))*nx/2;
    result[3] = -1.0*epsilon_*nx*nx+(a1_-fabs(a1_))*nx/2;
    result[2] = -1.0*epsilon_*ny*ny+(-a2_-fabs(a2_))*ny/2;
    result[4] = -1.0*epsilon_*ny*ny+(a2_-fabs(a2_))*ny/2;
    return result;
}

PositionArray ConDiff2D1::getJx(
	const Position,
	const Index,
	const Index ) const
{
    return jx_;
}

PositionArray ConDiff2D1::getJy(
	const Position,
	const Index,
	const Index ) const
{
    return jy_;
}

Index ConDiff2D1::size() const
{
    return 1;
}

bool ConDiff2D1::isConstant() const
{
    return true;
}

}
