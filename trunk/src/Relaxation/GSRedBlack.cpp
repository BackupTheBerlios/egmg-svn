/** \file GSRedBlack.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementaion of the class GSRedBlack.
 * \see GSRedBlack.h
 */
#include "GSRedBlack.h"

namespace mg
{

void GSRedBlack::relax(
    NumericArray& u,
    const NumericArray& f,
    const Stencil& stencil,
    const Index nx,
    const Index ny) const
{
    Precision factor = 1.0;
	if (stencil.size() < 2)
	{
		//first do the red points
		for (Index sy=1;sy<ny;sy+=2)
			for (Index sx=1;sx<nx;sx+=2)
			{
				factor = 1.0/stencil.getCenter(C,sx,sy,nx,ny);
				u[sy*(nx+1)+sx]+=factor*(f[sy*(nx+1)+sx]
							-stencil.apply(u,C,sx,sy,nx,ny));
			}
		for (Index sy=2;sy<(ny-1);sy+=2)
			for (Index sx=2;sx<(nx-1);sx+=2)
			{
				factor = 1.0/stencil.getCenter(C,sx,sy,nx,ny);
				u[sy*(nx+1)+sx]+=factor*(f[sy*(nx+1)+sx]
							-stencil.apply(u,C,sx,sy,nx,ny));
			}
		//do black points
		for (Index sy=1;sy<ny;sy+=2)
			for (Index sx=2;sx<(nx-1);sx+=2)
			{
				factor = 1.0/stencil.getCenter(C,sx,sy,nx,ny);
				u[sy*(nx+1)+sx]+=factor*(f[sy*(nx+1)+sx]
							-stencil.apply(u,C,sx,sy,nx,ny));
			}
		for (Index sy=2;sy<(ny-1);sy+=2)
			for (Index sx=1;sx<nx;sx+=2)
			{
				factor = 1.0/stencil.getCenter(C,sx,sy,nx,ny);
				u[sy*(nx+1)+sx]+=factor*(f[sy*(nx+1)+sx]
							-stencil.apply(u,C,sx,sy,nx,ny));
			}
	}
	else
	{
		/**
		 * \todo in case of large stencils a four colour RB is needed
		 */
		//first do the red points
		factor = 1.0/stencil.getCenter(SW,1,1,nx,ny);
		u[1*(nx+1)+1]+=factor*(f[1*(nx+1)+1]
					-stencil.apply(u,SW,1,1,nx,ny));
		for (Index sx=3;sx<(nx-1);sx+=2)
		{
			factor = 1.0/stencil.getCenter(S,sx,1,nx,ny);
			u[1*(nx+1)+sx]+=factor*(f[1*(nx+1)+sx]
					-stencil.apply(u,S,sx,1,nx,ny));
		}
		factor = 1.0/stencil.getCenter(SE,(nx-1),1,nx,ny);
		u[1*(nx+1)+(nx-1)]+=factor*(f[1*(nx+1)+(nx-1)]
					-stencil.apply(u,SE,(nx-1),1,nx,ny));
		for (Index sy=3;sy<(ny-1);sy+=2)
		{
			factor = 1.0/stencil.getCenter(W,1,sy,nx,ny);
			u[sy*(nx+1)+1]+=factor*(f[sy*(nx+1)+1]
					-stencil.apply(u,W,1,sy,nx,ny));
			for (Index sx=3;sx<(nx-1);sx+=2)
			{
				factor = 1.0/stencil.getCenter(C,sx,sy,nx,ny);
				u[sy*(nx+1)+sx]+=factor*(f[sy*(nx+1)+sx]
						-stencil.apply(u,C,sx,sy,nx,ny));
			}
		    factor = 1.0/stencil.getCenter(E,(nx-1),sy,nx,ny);
			u[sy*(nx+1)+(nx-1)]+=factor*(f[sy*(nx+1)+(nx-1)]
					-stencil.apply(u,E,(nx-1),sy,nx,ny));    
		}
		//the missing red points in the center
		for (Index sy=2;sy<(ny-1);sy+=2)
		{
			for (Index sx=2;sx<(nx-1);sx+=2)
			{
				factor = 1.0/stencil.getCenter(C,sx,sy,nx,ny);
				u[sy*(nx+1)+sx]+=factor*(f[sy*(nx+1)+sx]
						-stencil.apply(u,C,sx,sy,nx,ny));
			}
		}
		factor = 1.0/stencil.getCenter(NW,1,(ny-1),nx,ny);
		u[(nx-1)*(nx+1)+1]+=factor*(f[(nx-1)*(nx+1)+1]
					-stencil.apply(u,NW,1,(ny-1),nx,ny));
		for (Index sx=3;sx<(nx-1);sx+=2)
		{
			factor = 1.0/stencil.getCenter(N,sx,(nx-1),nx,ny);
			u[(nx-1)*(nx+1)+sx]+=factor*(f[(nx-1)*(nx+1)+sx]
					-stencil.apply(u,N,sx,(ny-1),nx,ny));
		}
		factor = 1.0/stencil.getCenter(NE,(nx-1),(nx-1),nx,ny);
		u[(nx-1)*(nx+1)+(nx-1)]+=factor*(f[(nx-1)*(nx+1)+(nx-1)]
				-stencil.apply(u,NE,(nx-1),(ny-1),nx,ny));
		
		//do black points
		for (Index sx=2;sx<(nx-1);sx+=2)
		{
			 factor = 1.0/stencil.getCenter(S,sx,1,nx,ny);
				u[1*(nx+1)+sx]+=factor*(f[1*(nx+1)+sx]
					-stencil.apply(u,S,sx,1,nx,ny));
		}
		for (Index sy=2;sy<(ny-1);sy+=2)
		{
			factor = 1.0/stencil.getCenter(W,1,sy,nx,ny);
			u[sy*(nx+1)+1]+=factor*(f[sy*(nx+1)+1]
					-stencil.apply(u,W,1,sy,nx,ny));	
			for (Index sx=3;sx<(nx-1);sx+=2)
			{
				factor = 1.0/stencil.getCenter(C,sx,sy,nx,ny);
				u[sy*(nx+1)+sx]+=factor*(f[sy*(nx+1)+sx]
						-stencil.apply(u,C,sx,sy,nx,ny));
			}
			factor = 1.0/stencil.getCenter(E,(nx-1),sy,nx,ny);
			u[sy*(nx+1)+(nx-1)]+=factor*(f[sy*(nx+1)+(nx-1)]
					-stencil.apply(u,E,(nx-1),sy,nx,ny));
		}
		for (Index sy=3;sy<(ny-1);sy+=2)
		{
			for (Index sx=2;sx<(nx-1);sx+=2)
			{
				factor = 1.0/stencil.getCenter(C,sx,sy,nx,ny);
				u[sy*(nx+1)+sx]+=factor*(f[sy*(nx+1)+sx]
						-stencil.apply(u,C,sx,sy,nx,ny));
			}
		}
		for (Index sx=2;sx<(nx-1);sx+=2)
		{
			factor = 1.0/stencil.getCenter(N,sx,(ny-1),nx,ny);
			u[(nx-1)*(nx+1)+sx]+=factor*(f[(nx-1)*(nx+1)+sx]
					-stencil.apply(u,N,sx,(ny-1),nx,ny));
		}
	}
}
}
