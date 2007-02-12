/** \file JumpingCoefficients.cpp
 * \author Kai Ruhnau, Jiri Kraus, Roman Wienands
 * \brief Contains the class implementation of JumpingCoefficients
 */

#include <cmath>

#include "JumpingCoefficients.h"
#include <iostream>
#include "../functions/printStencil.h"

namespace mg
{

PositionArray JumpingCoefficients::initJx_() const
{
    const int t[]={0,-1,0,1,0};
    return PositionArray(t,5);
}

PositionArray JumpingCoefficients::initJy_() const
{
    const int t[]={0,0,1,0,-1};
    return PositionArray(t,5);
}

Precision JumpingCoefficients::a_(
    const Precision x,
    const Precision y) const
{
    double h=1./128.0;
    switch ( mode_ )
    {
    case mA:
        if (x<0.5 && y<0.5)
            return 1;
        else if (x>=0.5 && y<0.5)
            return 1000;
        else if (x<0.5 && y>=0.5)
            return 10;
        else
            return 100;
    case mB:
        if (x<10*h || x>11*h)
            return 2;
        else
            return 1e-10;
    /*case mC:
        return -std::sin(Pi*sx/nx)*std::cos(Pi*sy/ny);
    case mD:
        return std::cos(15.0/180.0*Pi);*/
    default: 
        return 1;
    }
}

Precision JumpingCoefficients::b_(
    const Precision x,
    const Precision y) const
{
    double h=1./128.0;
    switch(mode_)
    {
    case mA:
        if (x<0.5 && y<0.5)
            return 1;
        else if (x>=0.5 && y<0.5)
            return 1000;
        else if (x<0.5 && y>=0.5)
            return 10;
        else
            return 100;
    case mB:
        if (x<10*h || x>11*h)
            return 2;
        else
            return 1e-10;
    /*case mC:
        return std::sin(Pi*sy/ny)*std::cos(Pi*sx/nx);
    case mD:
        return std::sin(15.0/180.0*Pi);*/
    default:
        return 1;
    }       
} 


JumpingCoefficients::JumpingCoefficients(Mode mode ) 
    : jx_( initJx_() ), jy_( initJy_() ), mode_( mode )
{}

JumpingCoefficients::~JumpingCoefficients() {}

Precision JumpingCoefficients::apply(
    const NumericArray& u,
    const Position,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny ) const
{
    return 
         ( nx*nx*(a_((sx-0.5)/nx,1.*sy/ny)+a_((sx+0.5)/nx,1.*sy/ny))
          +ny*ny*(b_(1.*sx/nx,(sy-0.5)/ny)+b_(1.*sx/nx,(sy+0.5)/ny))
         )*u[sy*(nx+1)+sx]

         -1.0*nx*nx * a_((sx-0.5)/nx,1.*sy/ny) * u[sy*(nx+1)+sx-1]
         -1.0*nx*nx * a_((sx+0.5)/nx,1.*sy/ny) * u[sy*(nx+1)+sx+1]
         -1.0*ny*ny * b_(1.*sx/nx,(sy-0.5)/ny) * u[(sy-1)*(nx+1)+sx]
         -1.0*ny*ny * b_(1.*sx/nx,(sy+0.5)/ny) * u[(sy+1)*(nx+1)+sx];

}

Precision JumpingCoefficients::getCenter(
   const Position,
   const Index sx,
   const Index sy,
   const Index nx,
   const Index ny ) const
{
    double value=
           nx*nx*(a_((sx-0.5)/nx,1.*sy/ny)+a_((sx+0.5)/nx,1.*sy/ny))
          +ny*ny*(b_(1.*sx/nx,(sy-0.5)/ny)+b_(1.*sx/nx,(sy+0.5)/ny));

    return value;
}

NumericArray JumpingCoefficients::getL(
    const Position,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny ) const
{
	NumericArray result( 0.0, 5 );
    result[C]=
           1.0*nx*nx*(a_((sx-0.5)/nx,1.*sy/ny)+a_((sx+0.5)/nx,1.*sy/ny))
          +1.0*ny*ny*(b_(1.*sx/nx,(sy-0.5)/ny)+b_(1.*sx/nx,(sy+0.5)/ny));

    result[W]=-1.0*nx*nx * a_((sx-0.5)/nx,1.*sy/ny);
    result[E]=-1.0*nx*nx * a_((sx+0.5)/nx,1.*sy/ny);
    result[S]=-1.0*ny*ny * b_(1.*sx/nx,(sy-0.5)/ny);
    result[N]=-1.0*ny*ny * b_(1.*sx/nx,(sy+0.5)/ny);
    return result; 
}

PositionArray JumpingCoefficients::getJx(
	const Position,
	const Index,
	const Index ) const
{
    return jx_;
}

PositionArray JumpingCoefficients::getJy(
	const Position,
	const Index,
	const Index ) const
{
    return jy_;
}

Index JumpingCoefficients::size() const
{
    return 1;
}

bool JumpingCoefficients::isConstant() const
{
    return false;
}

}
