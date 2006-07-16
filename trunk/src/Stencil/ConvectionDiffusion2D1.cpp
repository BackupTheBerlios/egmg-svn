/** \file ConvectionDiffusion2D1.cpp
 * \author Andre Oeckerath
 * \brief Contains the class implementation of ConvectionDiffusion2D1
 */

#include <cmath>

#include "ConvectionDiffusion2D1.h"

namespace mg
{

PositionArray ConvectionDiffusion2D1::initJx_() const
{
    const int t[]={0,-1,0,1,0};
    return PositionArray(t,5);
}

PositionArray ConvectionDiffusion2D1::initJy_() const
{
    const int t[]={0,0,1,0,-1};
    return PositionArray(t,5);
}

Precision ConvectionDiffusion2D1::a_(
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny) const
{
    switch ( mode_ )
    {
    case A:
        return (2.0*sy/ny-1)*(1-(sx*sx)/(nx*nx));
		break;
    case B:
        return 4.0*sx/nx*(sx/nx-1)*(1-2.0*sy/ny);
		break;
    default: 
        return 1;
    }
}

Precision ConvectionDiffusion2D1::b_(
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny) const
{
    switch(mode_)
    {
    case 0:
        return (2.0*sx*sy/(nx*ny))*(sy/ny-1);
    case 1:
        return -4.0*sy/ny*(sy/ny-1)*(1-2.0*sx/nx);
    default:
        return 1;
    }       
} 


ConvectionDiffusion2D1::ConvectionDiffusion2D1(Precision epsilon , Mode mode ) 
    : jx_( initJx_() ), jy_( initJy_() ),
	  epsilon_( epsilon ), mode_( mode )
{}

ConvectionDiffusion2D1::~ConvectionDiffusion2D1() {}

Precision ConvectionDiffusion2D1::apply(
    const NumericArray& u,
    const Position,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny ) const
{
    return 
         (2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny
          +fabs(a_(sx,sy,nx,ny))*nx+fabs(b_(sx,sy,nx,ny))*ny
         )*u[sy*(nx+1)+sx]
        +(-1.0*epsilon_*nx*nx
          +(-a_(sx,sy,nx,ny)-fabs(a_(sx,sy,nx,ny)))*nx/2
         )*u[sy*(nx+1)+sx-1]
        +(-1.0*epsilon_*nx*nx
          +(a_(sx,sy,nx,ny)-fabs(a_(sx,sy,nx,ny)))*nx/2
         )*u[sy*(nx+1)+sx+1]
        +(-1.0*epsilon_*ny*ny
          +(-b_(sx,sy,nx,ny)-fabs(b_(sx,sy,nx,ny)))*ny/2
         )*u[(sy-1)*(nx+1)+sx]
        +(-1.0*epsilon_*ny*ny
          +(b_(sx,sy,nx,ny)-fabs(b_(sx,sy,nx,ny)))*ny/2
         )*u[(sy+1)*(nx+1)+sx];
}

Precision ConvectionDiffusion2D1::getCenter(
   const Position,
   const Index sx,
   const Index sy,
   const Index nx,
   const Index ny ) const
{
    return 
         2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny
        +fabs( a_(sx,sy,nx,ny) )*nx+fabs( b_(sx,sy,nx,ny) )*ny;
}

NumericArray ConvectionDiffusion2D1::getL(
    const Position,
    const Index sx,
    const Index sy,
    const Index nx,
    const Index ny ) const
{
	NumericArray result( 0.0, 5 );
    result[0]=
         2.0*epsilon_*nx*nx+2.0*epsilon_*ny*ny
        +fabs(a_(sx,sy,nx,ny))*nx+fabs(b_(sx,sy,nx,ny))*ny;
    result[1]=
         -1.0*epsilon_*nx*nx
        +(-a_(sx,sy,nx,ny)-fabs(a_(sx,sy,nx,ny)))*nx/2;
    result[3]=
         -1.0*epsilon_*nx*nx
        +(a_(sx,sy,nx,ny)-fabs(a_(sx,sy,nx,ny)))*nx/2;
    result[2]=
         -1.0*epsilon_*ny*ny
        +(-b_(sx,sy,nx,ny)-fabs(b_(sx,sy,nx,ny)))*ny/2;
    result[4]=
         -1.0*epsilon_*ny*ny
        +(b_(sx,sy,nx,ny)-fabs(b_(sx,sy,nx,ny)))*ny/2;
    return result; 
}

PositionArray ConvectionDiffusion2D1::getJx(
	const Position,
	const Index,
	const Index ) const
{
    return jx_;
}

PositionArray ConvectionDiffusion2D1::getJy(
	const Position,
	const Index,
	const Index ) const
{
    return jy_;
}

Index ConvectionDiffusion2D1::size() const
{
    return 1;
}

bool ConvectionDiffusion2D1::isConstant() const
{
    return false;
}

}
