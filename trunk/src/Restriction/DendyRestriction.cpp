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
	//Dendy's restriction is defined here for max. compact 9-point stencils
	if (stencil.size() > 1)
		throw std::domain_error("u");

    const Index nxNew=nx/2;
    const Index nyNew=ny/2;
    NumericArray result(0.0,(nxNew+1)*(nyNew+1));
    NumericArray stencilL=NumericArray(9);
	Precision weightN; Precision weightS; Precision weightW; Precision weightE;
	Precision weightNW; Precision weightNE; Precision weightSW; Precision weightSE;

    //do injection on the boarders
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
			stencilL = stencil.getLinSize(C,2*sx,2*sy+1,nx,ny,1);
			weightN = -(stencilL[N] + stencilL[NW] + stencilL[NE]) / (stencilL[C] + stencilL[W] + stencilL[E]);
			stencilL = stencil.getLinSize(C,2*sx,2*sy-1,nx,ny,1);
			weightS = -(stencilL[S] + stencilL[SW] + stencilL[SE]) / (stencilL[C] + stencilL[W] + stencilL[E]);
			stencilL = stencil.getLinSize(C,2*sx+1,2*sy,nx,ny,1);
			weightE = -(stencilL[E] + stencilL[SE] + stencilL[NE]) / (stencilL[C] + stencilL[N] + stencilL[S]);
			stencilL = stencil.getLinSize(C,2*sx-1,2*sy,nx,ny,1);
			weightW = -(stencilL[W] + stencilL[SW] + stencilL[NW]) / (stencilL[C] + stencilL[N] + stencilL[S]);

			stencilL = stencil.getLinSize(C,2*sx-1,2*sy+1,nx,ny,1);		
			weightNW = -(stencilL[NW] + weightN * stencilL[W] + weightW * stencilL[N]) / stencilL[C];
			stencilL = stencil.getLinSize(C,2*sx+1,2*sy+1,nx,ny,1);			
			weightNE = -(stencilL[NE] + weightN * stencilL[E] + weightE * stencilL[N]) / stencilL[C];
			stencilL = stencil.getLinSize(C,2*sx-1,2*sy-1,nx,ny,1);			
			weightSW = -(stencilL[SW] + weightS * stencilL[W] + weightW * stencilL[S]) / stencilL[C];
			stencilL = stencil.getLinSize(C,2*sx+1,2*sy-1,nx,ny,1);			
			weightSE = -(stencilL[SE] + weightS * stencilL[E] + weightE * stencilL[S]) / stencilL[C];

            result[sy*(nxNew+1)+sx]=weight_ * (u[2*sy*(nx+1)+2*sx] 
				+ weightN*u[(2*sy+1)*(nx+1)+2*sx] + weightS*u[(2*sy-1)*(nx+1)+2*sx]
				+ weightW*u[2*sy*(nx+1)+2*sx-1] + weightE*u[2*sy*(nx+1)+2*sx+1]
				+ weightNW*u[(2*sy+1)*(nx+1)+2*sx-1] + weightNE*u[(2*sy+1)*(nx+1)+2*sx+1]
				+ weightSW*u[(2*sy-1)*(nx+1)+2*sx-1] + weightSE*u[(2*sy-1)*(nx+1)+2*sx+1]) / 4.0;
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
	if ( ( sx <= 1 ) || sx >= ( nx - 1 )  || ( sy <= 1 ) || ( sy >= ( ny - 1) ) )
	//if (sx==0 || sx==nx || sy==0 || sy==ny)
	{
		const Precision t[] = {
            1.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0};
		t_ = NumericArray(t,9);
	}
	else
	{
		NumericArray stencilL=NumericArray(9);
		t_.resize(9);
		stencilL = stencil.getLinSize(C,sx,sy+1,nx,ny,1);
		t_[N] = -(stencilL[N] + stencilL[NW] + stencilL[NE]) / (stencilL[C] + stencilL[W] + stencilL[E]);
		stencilL = stencil.getLinSize(C,sx,sy-1,nx,ny,1);
		t_[S] = -(stencilL[S] + stencilL[SW] + stencilL[SE]) / (stencilL[C] + stencilL[W] + stencilL[E]);
		stencilL = stencil.getLinSize(C,sx+1,sy,nx,ny,1);
		t_[E] = -(stencilL[E] + stencilL[SE] + stencilL[NE]) / (stencilL[C] + stencilL[N] + stencilL[S]);
		stencilL = stencil.getLinSize(C,sx-1,sy,nx,ny,1);
		t_[W] = -(stencilL[W] + stencilL[SW] + stencilL[NW]) / (stencilL[C] + stencilL[N] + stencilL[S]);

		stencilL = stencil.getLinSize(C,sx-1,sy+1,nx,ny,1);			
		t_[NW] = -(stencilL[NW] + t_[N] * stencilL[W] + t_[W] * stencilL[N]) / stencilL[C];
		stencilL = stencil.getLinSize(C,sx+1,sy+1,nx,ny,1);			
		t_[NE] = -(stencilL[NE] + t_[N] * stencilL[E] + t_[E] * stencilL[N]) / stencilL[C];
		stencilL = stencil.getLinSize(C,sx-1,sy-1,nx,ny,1);			
		t_[SW] = -(stencilL[SW] + t_[S] * stencilL[W] + t_[W] * stencilL[S]) / stencilL[C];
		stencilL = stencil.getLinSize(C,sx+1,sy-1,nx,ny,1);			
		t_[SE] = -(stencilL[SE] + t_[S] * stencilL[E] + t_[E] * stencilL[S]) / stencilL[C];

		t_[C]=4.0;

	}
	return t_;          
}
}

