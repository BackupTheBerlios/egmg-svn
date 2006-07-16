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
    PositionArray jx=stencil.getJx(C,nxNew,nyNew);
    PositionArray jy=stencil.getJy(C,nxNew,nyNew);
    NumericArray stencilL;
    std::valarray<Index> position( static_cast<Index>( 0 ), 9);

	for (Index j=0; j<jx.size(); ++j)
	{
		if (jy[j] == 1)		// north
		{
			if (jx[j] == -1)
				position[NW] = j;
			else if (jx[j] == 0)
				position[N] = j;
			else if (jx[j] == 1)
				position[NE] = j;
		}
		if (jy[j] == 0)		// center
		{
			if (jx[j] == -1)
				position[W] = j;
			else if (jx[j] == 0)
				position[C] = j;
			else if (jx[j] == 1)
				position[E] = j;
		}
		if (jy[j] == -1)	// south
		{
			if (jx[j] == -1)
				position[SW] = j;
			else if (jx[j] == 0)
				position[S] = j;
			else if (jx[j] == 1)
				position[SE] = j;
		}
	}

    Precision scale=0;
    Precision weight1=0;
    Precision weight2=0;
    Precision erg=0;

	// linear interpolation on the boarders
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
    
    //"interpolation" of coarse grid points
    for (Index sy=1; sy<=ny-1; sy++)
        for (Index sx=1; sx<=nx-1; sx++)
            result[2*sy*(nxNew+1)+2*sx]=u[sy*(nx+1)+sx];
    
	stencilL.resize( stencil.getL(C,2,1,nxNew,nyNew).size() );
    //interpolation of fine grid points on coarse grid lines
    for (Index sy=2; sy<=nyNew-2; sy+=2)
        for (Index sx=1; sx<=nxNew-1; sx+=2)
        {
            stencilL=stencil.getL(C,sx,sy,nxNew,nyNew);
            scale=-stencilL[0];
            weight1=0;
            weight2=0;
            if (position[S]!=0)
                scale-=stencilL[position[S]];
            if (position[N]!=0)
                scale-=stencilL[position[N]];
            if (position[W]!=0)
                weight1+=stencilL[position[W]];
            if (position[SW]!=0)
                weight1+=stencilL[position[SW]];
            if (position[NW]!=0)
                weight1+=stencilL[position[NW]];
            if (position[E]!=0)
                weight2+=stencilL[position[E]];
            if (position[SE]!=0)
                weight2+=stencilL[position[SE]];
            if (position[NE]!=0)
                weight2+=stencilL[position[NE]];
            result[sy*(nxNew+1)+sx]=(weight1*result[sy*(nxNew+1)+sx-1]+weight2*result[sy*(nxNew+1)+sx+1])/scale;
        }
    
    //interpolation of fine grid points on fine grid lines and coarse grid columns
    for (Index sy=1; sy<=nyNew-1; sy+=2)
        for (Index sx=2; sx<=nxNew-2; sx+=2)
        {
            stencilL=stencil.getL(C,sx,sy,nxNew,nyNew);
            scale=-stencilL[0];
            weight1=0;
            weight2=0;
            if (position[W]!=0)
                scale-=stencilL[position[W]];
            if (position[E]!=0)
                scale-=stencilL[position[E]];
            if (position[S]!=0)
                weight1+=stencilL[position[S]];
            if (position[SW]!=0)
                weight1+=stencilL[position[SW]];
            if (position[SE]!=0)
                weight1+=stencilL[position[SE]];
            if (position[N]!=0)
                weight2+=stencilL[position[N]];
            if (position[NW]!=0)
                weight2+=stencilL[position[NW]];
            if (position[NE]!=0)
                weight2+=stencilL[position[NE]];
            result[sy*(nxNew+1)+sx]=(weight1*result[(sy-1)*(nxNew+1)+sx]+weight2*result[(sy+1)*(nxNew+1)+sx])/scale;
        }
            
    //interpolation of fine grid points on fine grid lines and fine grid columns
    for (Index sy=1; sy<=nyNew-1; sy+=2)
        for (Index sx=1; sx<=nxNew-1; sx+=2)
        {
            stencilL=stencil.getL(C,sx,sy,nxNew,nyNew);
            erg=0;
            scale=-stencilL[0];
            if (position[W]!=0)
                erg+=stencilL[position[W]]*result[sy*(nxNew+1)+sx-1];
            if (position[E]!=0)
                erg+=stencilL[position[E]]*result[sy*(nxNew+1)+sx+1];
            if (position[S]!=0)
                erg+=stencilL[position[S]]*result[(sy-1)*(nxNew+1)+sx];
            if (position[SW]!=0)
                erg+=stencilL[position[SW]]*result[(sy-1)*(nxNew+1)+sx-1];
            if (position[SE]!=0)
                erg+=stencilL[position[SE]]*result[(sy-1)*(nxNew+1)+sx+1];
            if (position[N]!=0)
                erg+=stencilL[position[N]]*result[(sy+1)*(nxNew+1)+sx];
            if (position[NW]!=0)
                erg+=stencilL[position[NW]]*result[(sy+1)*(nxNew+1)+sx-1];
            if (position[NE]!=0)
                erg+=stencilL[position[NE]]*result[(sy+1)*(nxNew+1)+sx+1];
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
		PositionArray jx=stencil.getJx(C,nx,ny);
		PositionArray jy=stencil.getJy(C,nx,ny);
		std::valarray<Index> position((Index)0, 9);

		for (Index j=0; j<jx.size(); ++j)
		{
			if (jy[j] == 1)		// north
			{
				if (jx[j] == -1)
					position[NW] = j;
				else if (jx[j] == 0)
					position[N] = j;
				else if (jx[j] == 1)
					position[NE] = j;
			}
			if (jy[j] == 0)		// center
			{
				if (jx[j] == -1)
					position[W] = j;
				else if (jx[j] == 0)
					position[C] = j;
				else if (jx[j] == 1)
					position[E] = j;
			}
			if (jy[j] == -1)	// south
			{
				if (jx[j] == -1)
					position[SW] = j;
				else if (jx[j] == 0)
					position[S] = j;
				else if (jx[j] == 1)
					position[SE] = j;
			}
		}

		Precision scale=0;
		Precision weight1=0;
		Precision weight2=0;
		Precision erg=0;

		// C
		t_[0]=1.0;

		// W
		NumericArray stencilL=stencil.getL(C,sx-1,sy,nx,ny);
		scale=-stencilL[0];
		weight2=0;
		if (position[S]!=0)
			scale-=stencilL[position[S]];
		if (position[N]!=0)
			scale-=stencilL[position[N]];
		if (position[E]!=0)
			weight2+=stencilL[position[E]];
		if (position[SE]!=0)
			weight2+=stencilL[position[SE]];
		if (position[NE]!=0)
			weight2+=stencilL[position[NE]];
		t_[1]=weight2/scale; 

		// N
		stencilL=stencil.getL(C,sx,sy+1,nx,ny);
		scale=-stencilL[0];
		weight1=0;
		if (position[W]!=0)
			scale-=stencilL[position[W]];
		if (position[E]!=0)
			scale-=stencilL[position[E]];
		if (position[S]!=0)
			weight1+=stencilL[position[S]];
		if (position[SW]!=0)
			weight1+=stencilL[position[SW]];
		if (position[SE]!=0)
			weight1+=stencilL[position[SE]];
		t_[2]=weight1/scale;

		// E
		stencilL=stencil.getL(C,sx+1,sy,nx,ny);
		scale=-stencilL[0];
		weight1=0;
		if (position[S]!=0)
			scale-=stencilL[position[S]];
		if (position[N]!=0)
			scale-=stencilL[position[N]];
		if (position[W]!=0)
			weight1+=stencilL[position[W]];
		if (position[SW]!=0)
			weight1+=stencilL[position[SW]];
		if (position[NW]!=0)
			weight1+=stencilL[position[NW]];
		t_[3]=weight1/scale;      

		// S
		stencilL=stencil.getL(C,sx,sy-1,nx,ny);
		scale=-stencilL[0];
		weight2=0;
		if (position[W]!=0)
			scale-=stencilL[position[W]];
		if (position[E]!=0)
			scale-=stencilL[position[E]];
		if (position[N]!=0)
			weight2+=stencilL[position[N]];
		if (position[NW]!=0)
			weight2+=stencilL[position[NW]];
		if (position[NE]!=0)
			weight2+=stencilL[position[NE]];
		t_[4]=weight2/scale;

		// NW
		stencilL=stencil.getL(C,sx-1,sy+1,nx,ny);
		scale=-stencilL[0];
		erg=0;
		if (position[E]!=0)
			erg+=stencilL[position[E]]*t_[2];
		if (position[S]!=0)
			erg+=stencilL[position[S]]*t_[1];
		if (position[SE]!=0)
			erg+=stencilL[position[SE]];
		t_[5]=erg/scale;

		// NE
		stencilL=stencil.getL(C,sx+1,sy+1,nx,ny);
		scale=-stencilL[0];
		erg=0;
		if (position[W]!=0)
			erg+=stencilL[position[W]]*t_[2];
		if (position[S]!=0)
			erg+=stencilL[position[S]]*t_[3];
		if (position[SW]!=0)
			erg+=stencilL[position[SW]];
		t_[6]=erg/scale;

		// SE
		stencilL=stencil.getL(C,sx+1,sy-1,nx,ny);
		scale=-stencilL[0];
		erg=0;
		if (position[W]!=0)
			erg+=stencilL[position[W]]*t_[4];
		if (position[N]!=0)
			erg+=stencilL[position[N]]*t_[3];
		if (position[NW]!=0)
			erg+=stencilL[position[NW]];
		t_[7]=erg/scale;

		// SW
		stencilL=stencil.getL(C,sx-1,sy-1,nx,ny);
		scale=-stencilL[0];
		erg=0;
		if (position[E]!=0)
			erg+=stencilL[position[E]]*t_[4];
		if (position[N]!=0)
			erg+=stencilL[position[N]]*t_[1];
		if (position[NE]!=0)
			erg+=stencilL[position[NE]];
		t_[8]=erg/scale;      

		return t_;          
    }
    const Precision t[] = {
            1.0,   1.0/2, 1.0/2,
            1.0/2, 1.0/2, 1.0/4,
            1.0/4, 1.0/4, 1.0/4};
    t_ = NumericArray(t,9);
    return t_;
}

}
