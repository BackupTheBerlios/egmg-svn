/** \file DendyInterpolation.cpp
 * \author Benedikt Engbroks
 * \brief
 Contains the implementation of the class DendyInterpolation
 * \see DendyInterpolation.h
 */
#include "DendyInterpolation.h"

namespace mg
{
// This function only works for max. compact 9-point stencils
NumericArray DendyInterpolation::prolongate(
    const NumericArray& u,
    const Stencil& stencil,
    const Index nx,
    const Index ny) const
{
    const Index nxNew=2*nx;
    const Index nyNew=2*ny;
    NumericArray result((nxNew+1)*(nyNew+1));
    NumericArray stencilL;

    Precision scale=0;
    Precision weight1=0;
    Precision weight2=0;
    Precision erg=0;

	// linear interpolation on the borders
	for (Index sy=0; sy<=ny; sy++)
	{
		result[2*sy*(nxNew+1)]=u[sy*(nx+1)];
        result[2*sy*(nxNew+1)+2*nx]=u[sy*(nx+1)+nx];
    }

	for (Index sx=0; sx<=nx; sx++)
	{
		result[2*sx]=u[sx];
		result[2*ny*(nxNew+1)+2*sx]=u[ny*(nx+1)+sx];
	}

	for (Index sy=0; sy<=ny-1; sy++)
	{
		result[(2*sy+1)*(nxNew+1)] = 0.5 * (result[2*sy*(nxNew+1)] + result[2*(sy+1)*(nxNew+1)]);
		result[(2*sy+1)*(nxNew+1)+2*nx] = 0.5 * (result[2*sy*(nxNew+1)+2*nx] + result[2*(sy+1)*(nxNew+1)+2*nx]);
	}

	for (Index sx=0; sx<=nx-1; sx++)
	{
		result[2*sx+1] = 0.5 * (result[2*sx] + result[2*sx+2]);
		result[2*ny*(nxNew+1)+2*sx+1] = 0.5 * (result[2*ny*(nxNew+1)+2*sx] + result[2*ny*(nxNew+1)+2*sx+2]);
	}
    
    //copy coarse grid points
    for (Index sy=1; sy<=ny-1; sy++)
        for (Index sx=1; sx<=nx-1; sx++)
            result[2*sy*(nxNew+1)+2*sx]=u[sy*(nx+1)+sx];
    
	stencilL.resize( stencil.getLInSize(C,2,1,nxNew,nyNew, 1).size());
    //interpolation of fine grid points on coarse grid lines
    for (Index sy=2; sy<=nyNew-2; sy+=2)
        for (Index sx=1; sx<=nxNew-1; sx+=2)
        {
            stencilL=stencil.getLInSize(C,sx,sy,nxNew,nyNew, 1);

            scale=0;
            scale-=stencilL[C];
            scale-=stencilL[S];
            scale-=stencilL[N];

            weight1=0;
            weight1+=stencilL[W];
            weight1+=stencilL[SW];
            weight1+=stencilL[NW];

            weight2=0;
            weight2+=stencilL[E];
            weight2+=stencilL[SE];
            weight2+=stencilL[NE];

            result[sy*(nxNew+1)+sx]=
                    ( weight1*result[sy*(nxNew+1)+sx-1]
                     +weight2*result[sy*(nxNew+1)+sx+1]
                    )/scale;
        }
    
    //interpolation of fine grid points on fine grid lines and coarse grid columns
    for (Index sy=1; sy<=nyNew-1; sy+=2)
        for (Index sx=2; sx<=nxNew-2; sx+=2)
        {
            stencilL=stencil.getLInSize(C,sx,sy,nxNew,nyNew,1);

            scale=0;
            scale-=stencilL[C];
            scale-=stencilL[W];
            scale-=stencilL[E];

            weight1=0;
            weight1+=stencilL[S];
            weight1+=stencilL[SW];
            weight1+=stencilL[SE];

            weight2=0;
            weight2+=stencilL[N];
            weight2+=stencilL[NW];
            weight2+=stencilL[NE];

            result[sy*(nxNew+1)+sx]=
                    ( weight1*result[(sy-1)*(nxNew+1)+sx]
                     +weight2*result[(sy+1)*(nxNew+1)+sx]
                    )/scale;
        }
            
    //interpolation of fine grid points on fine grid lines and fine grid columns
    for (Index sy=1; sy<=nyNew-1; sy+=2)
        for (Index sx=1; sx<=nxNew-1; sx+=2)
        {
            stencilL=stencil.getLInSize(C,sx,sy,nxNew,nyNew,1);
            erg=0;
            scale=-stencilL[0];

            erg+=stencilL[W] *result[ sy   *(nxNew+1)+sx-1];
            erg+=stencilL[E] *result[ sy   *(nxNew+1)+sx+1];
            erg+=stencilL[S] *result[(sy-1)*(nxNew+1)+sx  ];
            erg+=stencilL[SW]*result[(sy-1)*(nxNew+1)+sx-1];
            erg+=stencilL[SE]*result[(sy-1)*(nxNew+1)+sx+1];
            erg+=stencilL[N] *result[(sy+1)*(nxNew+1)+sx  ];
            erg+=stencilL[NW]*result[(sy+1)*(nxNew+1)+sx-1];
            erg+=stencilL[NE]*result[(sy+1)*(nxNew+1)+sx+1];

            result[sy*(nxNew+1)+sx]=erg/scale;
        }

    return result;
}

const NumericArray& DendyInterpolation::getI(
	const Position, 
    const Index sx,
    const Index sy, 
    const Index nx,
    const Index ny,
    const Stencil& stencil) const
{
	if ( ! ( ( sx <= 1 ) || sx >= ( nx - 1 )  || ( sy <= 1 ) || ( sy >= ( ny - 1) ) ) )
	{
		Precision scale=0;
		Precision weight1=0;
		Precision weight2=0;
		Precision erg=0;

		t_[C]=1.0;

		NumericArray stencilL=stencil.getLInSize(C,sx-1,sy,nx,ny,1);

        scale=0;
		scale-=stencilL[C];
		scale-=stencilL[S];
		scale-=stencilL[N];

		weight2=0;
		weight2+=stencilL[E];
		weight2+=stencilL[SE];
		weight2+=stencilL[NE];

		t_[W]=weight2/scale; 

		stencilL=stencil.getLInSize(C,sx,sy+1,nx,ny,1);

        scale=0;
		scale-=stencilL[C];
		scale-=stencilL[W];
		scale-=stencilL[E];

		weight1=0;
		weight1+=stencilL[S];
		weight1+=stencilL[SW];
		weight1+=stencilL[SE];

		t_[N]=weight1/scale;

		stencilL=stencil.getLInSize(C,sx+1,sy,nx,ny,1);

        scale=0;
		scale-=stencilL[C];
		scale-=stencilL[S];
		scale-=stencilL[N];

		weight1=0;
		weight1+=stencilL[W];
		weight1+=stencilL[SW];
		weight1+=stencilL[NW];

		t_[E]=weight1/scale;      

		stencilL=stencil.getLInSize(C,sx,sy-1,nx,ny,1);

        scale=0;
		scale-=stencilL[C];
		scale-=stencilL[W];
		scale-=stencilL[E];

		weight2=0;
		weight2+=stencilL[N];
		weight2+=stencilL[NW];
		weight2+=stencilL[NE];

		t_[S]=weight2/scale;

		stencilL=stencil.getLInSize(C,sx-1,sy+1,nx,ny,1);
		scale=-stencilL[C];

		erg=0;
		erg+=stencilL[SE];
		erg+=stencilL[E]*t_[N];
		erg+=stencilL[S]*t_[W];

		t_[NW]=erg/scale;

		stencilL=stencil.getLInSize(C,sx+1,sy+1,nx,ny,1);
		scale=-stencilL[C];

		erg=0;
		erg+=stencilL[SW];
		erg+=stencilL[W]*t_[N];
		erg+=stencilL[S]*t_[E];

		t_[NE]=erg/scale;

		stencilL=stencil.getLInSize(C,sx+1,sy-1,nx,ny,1);
		scale=-stencilL[C];

		erg=0;
		erg+=stencilL[NW];
		erg+=stencilL[W]*t_[S];
		erg+=stencilL[N]*t_[E];

		t_[SE]=erg/scale;

		stencilL=stencil.getLInSize(C,sx-1,sy-1,nx,ny,1);
		scale=-stencilL[C];

		erg=0;
		erg+=stencilL[NE];
		erg+=stencilL[E]*t_[S];
		erg+=stencilL[N]*t_[W];

		t_[SW]=erg/scale;      
    }
    else
    {
        const Precision t[] = {
                1.0,   1.0/2, 1.0/2,
                1.0/2, 1.0/2, 1.0/4,
                1.0/4, 1.0/4, 1.0/4};
        t_ = NumericArray(t,9);
    }
    return t_;
}

}
