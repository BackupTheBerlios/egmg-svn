/** \file DendyRestriction.cpp
 * \author Benedikt Engbroks
 * \brief Contains the implementation of the class DendyRestriction.
 */

#include "DendyRestriction.h"
#include <stdexcept>

namespace mg
{
NumericArray DendyRestriction::restriction(
    const NumericArray& u,
    const Stencil& stencil,
    const Index nx,const Index ny) const
{
    //if it is not possible to do standart coarsening throw an exeption
    if ((nx%2 != 0) || (ny%2 != 0))
        throw std::domain_error("u");

    const Index nxNew=nx/2;
    const Index nyNew=ny/2;
    NumericArray result(0.0,(nxNew+1)*(nyNew+1));
    NumericArray stencilL=NumericArray(9);
    Precision weight[NamedPositions]={0};

    //do injection on the borders
    for (Index sy=0;sy<=nyNew;sy++)
    {
        result[sy*(nxNew+1)]=u[2*sy*(nx+1)];
        result[sy*(nxNew+1)+nxNew]=u[2*sy*(nx+1)+nx];
    }
    for (Index sx=0;sx<=nxNew;sx++)
    {
        result[sx]=u[2*sx];
        result[nyNew*(nxNew+1)+sx]=u[ny*(nx+1)+2*sx];
    }

    for (Index sy=1; sy<=nyNew-1; sy++)
        for (Index sx=1; sx<=nxNew-1; sx++)
		{
			stencilL = stencil.getLInSize(C,2*sx,2*sy+1,nx,ny,1);
			weight[N] = -(stencilL[N] + stencilL[NW] + stencilL[NE])
                       / (stencilL[C] + stencilL[W]  + stencilL[E]);
			stencilL = stencil.getLInSize(C,2*sx,2*sy-1,nx,ny,1);
			weight[S] = -(stencilL[S] + stencilL[SW] + stencilL[SE]) 
                       / (stencilL[C] + stencilL[W]  + stencilL[E]);
			stencilL = stencil.getLInSize(C,2*sx+1,2*sy,nx,ny,1);
			weight[E] = -(stencilL[E] + stencilL[SE] + stencilL[NE])
                       / (stencilL[C] + stencilL[N]  + stencilL[S]);
			stencilL = stencil.getLInSize(C,2*sx-1,2*sy,nx,ny,1);
			weight[W] = -(stencilL[W] + stencilL[SW] + stencilL[NW])
                       / (stencilL[C] + stencilL[N]  + stencilL[S]);

			stencilL = stencil.getLInSize(C,2*sx-1,2*sy+1,nx,ny,1);		
			weight[NW] = -(  stencilL[NW] 
                           + weight[N] * stencilL[W] 
                           + weight[W] * stencilL[N]) / stencilL[C];
			stencilL = stencil.getLInSize(C,2*sx+1,2*sy+1,nx,ny,1);			
			weight[NE] = -(  stencilL[NE] 
                           + weight[N] * stencilL[E]
                           + weight[E] * stencilL[N]) / stencilL[C];
			stencilL = stencil.getLInSize(C,2*sx-1,2*sy-1,nx,ny,1);			
			weight[SW] = -(  stencilL[SW] 
                           + weight[S] * stencilL[W] 
                           + weight[W] * stencilL[S]) / stencilL[C];
			stencilL = stencil.getLInSize(C,2*sx+1,2*sy-1,nx,ny,1);			
			weight[SE] = -(  stencilL[SE]
                           + weight[S] * stencilL[E] 
                           + weight[E] * stencilL[S]) / stencilL[C];

            result[sy*(nxNew+1)+sx]=weight_/4.0 * (
                             u[2*sy*(nx+1)+2*sx] 
				+ weight[N]* u[(2*sy+1)*(nx+1)+2*sx] 
                + weight[S]*u[(2*sy-1)*(nx+1)+2*sx]
				+ weight[W]* u[2*sy*(nx+1)+2*sx-1] 
                + weight[E]*u[2*sy*(nx+1)+2*sx+1]
				+ weight[NW]*u[(2*sy+1)*(nx+1)+2*sx-1] 
                + weight[NE]*u[(2*sy+1)*(nx+1)+2*sx+1]
				+ weight[SW]*u[(2*sy-1)*(nx+1)+2*sx-1]
                + weight[SE]*u[(2*sy-1)*(nx+1)+2*sx+1]);
		}
    
    return result;
}

const NumericArray& DendyRestriction::getI(
	const Position, 
    const Index sx,
    const Index sy, 
    const Index nx,
    const Index ny,
    const Stencil& stencil) const
{
	if (   sx <= 1 || sx >= ( nx - 1 ) 
        || sy <= 1 || sy >= ( ny - 1 ) )
	{
		const Precision t[] = {
            1.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0};
		t_.resize(9);
		t_ = NumericArray(t,9);
	}
	else
	{
		NumericArray stencilL=NumericArray(9);
		t_.resize(9);
		stencilL = stencil.getLInSize(C,sx,sy+1,nx,ny,1);
		t_[N] = -(stencilL[N] + stencilL[NW] + stencilL[NE]) 
                /(stencilL[C] + stencilL[W]  + stencilL[E] );
		stencilL = stencil.getLInSize(C,sx,sy-1,nx,ny,1);
		t_[S] = -(stencilL[S] + stencilL[SW] + stencilL[SE])
                /(stencilL[C] + stencilL[W]  + stencilL[E] );
		stencilL = stencil.getLInSize(C,sx+1,sy,nx,ny,1);
		t_[E] = -(stencilL[E] + stencilL[SE] + stencilL[NE])
                /(stencilL[C] + stencilL[N]  + stencilL[S] );
		stencilL = stencil.getLInSize(C,sx-1,sy,nx,ny,1);
		t_[W] = -(stencilL[W] + stencilL[SW] + stencilL[NW])
                /(stencilL[C] + stencilL[N]  + stencilL[S] );

		stencilL = stencil.getLInSize(C,sx-1,sy+1,nx,ny,1);			
		t_[NW] = -(stencilL[NW] + t_[N] * stencilL[W] + t_[W] * stencilL[N]) 
                 / stencilL[C];
		stencilL = stencil.getLInSize(C,sx+1,sy+1,nx,ny,1);			
		t_[NE] = -(stencilL[NE] + t_[N] * stencilL[E] + t_[E] * stencilL[N])
                 / stencilL[C];
		stencilL = stencil.getLInSize(C,sx-1,sy-1,nx,ny,1);			
		t_[SW] = -(stencilL[SW] + t_[S] * stencilL[W] + t_[W] * stencilL[S])
                 / stencilL[C];
		stencilL = stencil.getLInSize(C,sx+1,sy-1,nx,ny,1);			
		t_[SE] = -(stencilL[SE] + t_[S] * stencilL[E] + t_[E] * stencilL[S])
                 / stencilL[C];

		t_[C]=4.0;
	}

	return t_;          
}
}

