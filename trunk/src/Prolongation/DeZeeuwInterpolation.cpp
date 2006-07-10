/** \file DeZeeuwInterpolation.cpp
 * \author Benedikt Engbroks
 * \brief Contains the implementation of the class DeZeeuwInterpolation
 * \see DeZeeuwInterpolation.h
 * \todo This file needs cleanup
 */

#include "DeZeeuwInterpolation.h"
#include "../Stencil/Stencil.h"
#include <algorithm>

namespace mg
{
// This function only works for max. compact 9-point stencils
NumericArray DeZeeuwInterpolation::prolongate(
    const NumericArray& u,
    const Stencil& stencil,
    const Index nx,
    const Index ny) const
{
    register const Index nxNew=2*nx;
    register const Index nyNew=2*ny;
    NumericArray result((nxNew+1)*(nyNew+1));
    PositionArray jx=stencil.getJx(C,nx,ny);
    PositionArray jy=stencil.getJy(C,nx,ny);
    NumericArray stencilL=stencil.getL(C,0,0,nx,ny);
    std::valarray<Index> position((Index)0, 9);
    NumericArray ms(9);
    NumericArray mt(9);
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
    Precision symsum=0;
    Precision d_w=0;
    Precision d_e=0;
    Precision d_n=0;
    Precision d_s=0;
    Precision sigma1=0;
    Precision c_1=0;
    Precision w_w=0;
    Precision w_e=0;
    Precision sigma2=0;
    Precision c_2=0;
    Precision w_s=0;
    Precision w_n=0;

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
    
    //interpolation of fine grid points on coarse grid lines
    for (Index sy=2; sy<=nyNew-2; sy+=2)
        for (Index sx=1; sx<=nxNew-1; sx+=2)
        {
            stencilL=stencil.getL(C,sx,sy,nx,ny);
            symsum=0;
            
            // Divide the stencil defined by stencilL und position into a
            // symmetric and an antisymmetric part.

			ms[NE]=ms[SW]=0.5*(stencilL[position[SW]]+stencilL[position[NE]]);
			ms[N]=ms[S]=0.5*(stencilL[position[S]]+stencilL[position[N]]);
			ms[NW]=ms[SE]=0.5*(stencilL[position[SE]]+stencilL[position[NW]]);
			ms[E]=ms[W]=0.5*(stencilL[position[W]]+stencilL[position[E]]);
			ms[C]=stencilL[position[C]];
            symsum+=2*ms[NE]+2*ms[N]+2*ms[NW]+2*ms[E]+ms[C];
			mt[SW]=0.5*(stencilL[position[SW]]-stencilL[position[NE]]);
			mt[S]=0.5*(stencilL[position[S]]-stencilL[position[N]]);
			mt[SE]=0.5*(stencilL[position[SE]]-stencilL[position[NW]]);
			mt[W]=0.5*(stencilL[position[W]]-stencilL[position[E]]);
			mt[NE]=0.5*(stencilL[position[NE]]-stencilL[position[SW]]);
			mt[N]=0.5*(stencilL[position[N]]-stencilL[position[S]]);
			mt[NW]=0.5*(stencilL[position[NW]]-stencilL[position[SE]]);
			mt[E]=0.5*(stencilL[position[E]]-stencilL[position[E]]);
			mt[C]=0.0;

			d_w=std::max(std::fabs(ms[SW]+ms[W]+ms[NW]),
                         std::max(std::fabs(ms[SW]),
                                  std::fabs(ms[NW])));
            d_e=std::max(std::fabs(ms[SE]+ms[E]+ms[NE]),
                         std::max(std::fabs(ms[SE]),
                                  std::fabs(ms[NE])));
            d_n=std::max(std::fabs(ms[NW]+ms[N]+ms[NE]),
                         std::max(std::fabs(ms[NW]),
                                  std::fabs(ms[NE])));
            d_s=std::max(std::fabs(ms[SW]+ms[S]+ms[SE]),
                         std::max(std::fabs(ms[SW]),
                                  std::fabs(ms[SE])));
            sigma1=0.5*std::min(1.0, std::fabs(1-symsum/stencilL[0]));
            c_1=mt[SE]+mt[E]+mt[NE]-mt[SW]-mt[W]-mt[NW];
            w_w=sigma1*(1+(d_w-d_e)/(d_w+d_e)+c_1/(d_w+d_e+d_n+d_s));
            w_e=2*sigma1-w_w;
            weight1=std::min(2*sigma1, std::max(w_w, 0.0));
            weight2=std::min(2*sigma1, std::max(w_e, 0.0));
            result[sy*(nxNew+1)+sx]=
                 weight1*result[sy*(nxNew+1)+sx-1]
                +weight2*result[sy*(nxNew+1)+sx+1];
        }
    
    // interpolation of fine grid points on fine grid lines and coarse grid 
    // columns
    for (Index sy=1; sy<=nyNew-1; sy+=2)
        for (Index sx=2; sx<=nxNew-2; sx+=2)
        {
            stencilL=stencil.getL(C,sx,sy,nx,ny);
            symsum=0;
            
            // Divide the stencil defined by stencilL und position into a
            // symmetric and an antisymmetric part.
			ms[NE]=ms[SW]=0.5*(stencilL[position[SW]]+stencilL[position[NE]]);
			ms[N]=ms[S]=0.5*(stencilL[position[S]]+stencilL[position[N]]);
			ms[NW]=ms[SE]=0.5*(stencilL[position[SE]]+stencilL[position[NW]]);
			ms[E]=ms[W]=0.5*(stencilL[position[W]]+stencilL[position[E]]);
			ms[C]=stencilL[position[C]];
            symsum+=2*ms[NE]+2*ms[N]+2*ms[NW]+2*ms[E]+ms[C];
			mt[SW]=0.5*(stencilL[position[SW]]-stencilL[position[NE]]);
			mt[S]=0.5*(stencilL[position[S]]-stencilL[position[N]]);
			mt[SE]=0.5*(stencilL[position[SE]]-stencilL[position[NW]]);
			mt[W]=0.5*(stencilL[position[W]]-stencilL[position[E]]);
			mt[NE]=0.5*(stencilL[position[NE]]-stencilL[position[SW]]);
			mt[N]=0.5*(stencilL[position[N]]-stencilL[position[S]]);
			mt[NW]=0.5*(stencilL[position[NW]]-stencilL[position[SE]]);
			mt[E]=0.5*(stencilL[position[E]]-stencilL[position[E]]);
			mt[C]=0.0;

			d_w=std::max(std::fabs(ms[SW]+ms[W]+ms[NW]),
                         std::max(std::fabs(ms[SW]),
                                  std::fabs(ms[NW])));
            d_e=std::max(std::fabs(ms[SE]+ms[E]+ms[NE]),
                         std::max(std::fabs(ms[SE]),
                                  std::fabs(ms[NE])));
            d_n=std::max(std::fabs(ms[NW]+ms[N]+ms[NE]),
                         std::max(std::fabs(ms[NW]),
                                  std::fabs(ms[NE])));
            d_s=std::max(std::fabs(ms[SW]+ms[S]+ms[SE]),
                         std::max(std::fabs(ms[SW]),
                                  std::fabs(ms[SE])));
            sigma2=0.5*std::min(1.0, std::fabs(1-symsum/stencilL[0]));
            c_2=mt[NW]+mt[N]+mt[NE]-mt[SW]-mt[S]-mt[SE];
            w_n=sigma2*(1+(d_s-d_n)/(d_s+d_n)+c_2/(d_w+d_e+d_n+d_s));
            w_s=2*sigma2-w_n;
            weight1=std::min(2*sigma2, std::max(w_s, 0.0));
            weight2=std::min(2*sigma2, std::max(w_n, 0.0));
            result[sy*(nxNew+1)+sx]=
                 weight1*result[(sy-1)*(nxNew+1)+sx]
                +weight2*result[(sy+1)*(nxNew+1)+sx];
        }
            
    //interpolation of fine grid points on fine grid lines and fine grid columns
    for (Index sy=1; sy<=nyNew-1; sy+=2)
        for (Index sx=1; sx<=nxNew-1; sx+=2)
        {
            stencilL=stencil.getL(C,sx,sy,nx,ny);
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
const NumericArray& DeZeeuwInterpolation::getI(
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
		NumericArray ms(9);
		NumericArray mt(9);

		Precision scale=0;
		Precision weight1=0;
		Precision weight2=0;
		Precision erg=0;
		Precision symsum=0;
		Precision d_w=0;
		Precision d_e=0;
		Precision d_n=0;
		Precision d_s=0;
		Precision sigma1=0;
		Precision c_1=0;
		Precision w_w=0;
		Precision w_e=0;
		Precision sigma2=0;
		Precision c_2=0;
		Precision w_s=0;
		Precision w_n=0;

		// C
		t_[0]=1.0;

		// W
		NumericArray stencilL=stencil.getL(C,sx-1,sy,nx,ny);
		symsum=0;
	    
		// Divide the stencil defined by stencilL und position into a symmetric 
		// and an antisymmetric part.
		ms[NE]=ms[SW]=0.5*(stencilL[position[SW]]+stencilL[position[NE]]);
		ms[N]=ms[S]=0.5*(stencilL[position[S]]+stencilL[position[N]]);
		ms[NW]=ms[SE]=0.5*(stencilL[position[SE]]+stencilL[position[NW]]);
		ms[E]=ms[W]=0.5*(stencilL[position[W]]+stencilL[position[E]]);
		ms[C]=stencilL[position[C]];
		symsum+=2*ms[NE]+2*ms[N]+2*ms[NW]+2*ms[E]+ms[C];
		mt[SW]=0.5*(stencilL[position[SW]]-stencilL[position[NE]]);
		mt[S]=0.5*(stencilL[position[S]]-stencilL[position[N]]);
		mt[SE]=0.5*(stencilL[position[SE]]-stencilL[position[NW]]);
		mt[W]=0.5*(stencilL[position[W]]-stencilL[position[E]]);
		mt[NE]=0.5*(stencilL[position[NE]]-stencilL[position[SW]]);
		mt[N]=0.5*(stencilL[position[N]]-stencilL[position[S]]);
		mt[NW]=0.5*(stencilL[position[NW]]-stencilL[position[SE]]);
		mt[E]=0.5*(stencilL[position[E]]-stencilL[position[E]]);
		mt[C]=0.0;
		d_w=std::max(std::fabs(ms[SW]+ms[W]+ms[NW]),
					 std::max(std::fabs(ms[SW]),
							  std::fabs(ms[NW])));
		d_e=std::max(std::fabs(ms[SE]+ms[E]+ms[NE]),
					 std::max(std::fabs(ms[SE]),
							  std::fabs(ms[NE])));
		d_n=std::max(std::fabs(ms[NW]+ms[N]+ms[NE]),
					 std::max(std::fabs(ms[NW]),
							  std::fabs(ms[NE])));
		d_s=std::max(std::fabs(ms[SW]+ms[S]+ms[SE]),
					 std::max(std::fabs(ms[SW]),
							  std::fabs(ms[SE])));
		sigma1=0.5*std::min(1.0, std::fabs(1-symsum/stencilL[0]));
		c_1=mt[SE]+mt[E]+mt[NE]-mt[SW]-mt[W]-mt[NW];
		w_w=sigma1*(1+(d_w-d_e)/(d_w+d_e)+c_1/(d_w+d_e+d_n+d_s));
		w_e=2*sigma1-w_w;
		weight2=std::min(2*sigma1, std::max(w_e, 0.0));
		t_[1]=weight2; 

		// N
		stencilL=stencil.getL(C,sx,sy+1,nx,ny);
		symsum=0;
	    
		// Divide the stencil defined by stencilL und position into a symmetric 
		// and an antisymmetric part.
		ms[NE]=ms[SW]=0.5*(stencilL[position[SW]]+stencilL[position[NE]]);
		ms[N]=ms[S]=0.5*(stencilL[position[S]]+stencilL[position[N]]);
		ms[NW]=ms[SE]=0.5*(stencilL[position[SE]]+stencilL[position[NW]]);
		ms[E]=ms[W]=0.5*(stencilL[position[W]]+stencilL[position[E]]);
		ms[C]=stencilL[position[C]];
		symsum+=2*ms[NE]+2*ms[N]+2*ms[NW]+2*ms[E]+ms[C];
		mt[SW]=0.5*(stencilL[position[SW]]-stencilL[position[NE]]);
		mt[S]=0.5*(stencilL[position[S]]-stencilL[position[N]]);
		mt[SE]=0.5*(stencilL[position[SE]]-stencilL[position[NW]]);
		mt[W]=0.5*(stencilL[position[W]]-stencilL[position[E]]);
		mt[NE]=0.5*(stencilL[position[NE]]-stencilL[position[SW]]);
		mt[N]=0.5*(stencilL[position[N]]-stencilL[position[S]]);
		mt[NW]=0.5*(stencilL[position[NW]]-stencilL[position[SE]]);
		mt[E]=0.5*(stencilL[position[E]]-stencilL[position[E]]);
		mt[C]=0.0;
		d_w=std::max(std::fabs(ms[SW]+ms[W]+ms[NW]),
					 std::max(std::fabs(ms[SW]),
							  std::fabs(ms[NW])));
		d_e=std::max(std::fabs(ms[SE]+ms[E]+ms[NE]),
					 std::max(std::fabs(ms[SE]),
							  std::fabs(ms[NE])));
		d_n=std::max(std::fabs(ms[NW]+ms[N]+ms[NE]),
					 std::max(std::fabs(ms[NW]),
							  std::fabs(ms[NE])));
		d_s=std::max(std::fabs(ms[SW]+ms[S]+ms[SE]),
					 std::max(std::fabs(ms[SW]),
							  std::fabs(ms[SE])));
		sigma2=0.5*std::min(1.0, std::fabs(1-symsum/stencilL[0]));
		c_2=mt[NW]+mt[N]+mt[NE]-mt[SW]-mt[S]-mt[SE];
		w_n=sigma2*(1+(d_s-d_n)/(d_s+d_n)+c_2/(d_w+d_e+d_n+d_s));
		w_s=2*sigma2-w_n;
		weight1=std::min(2*sigma2, std::max(w_s, 0.0));
		weight2=std::min(2*sigma2, std::max(w_n, 0.0));
		t_[2]=weight1;

		// E
		stencilL=stencil.getL(C,sx+1,sy,nx,ny);
		symsum=0;
	    
		// Divide the stencil defined by stencilL und position into a symmetric 
		// and an antisymmetric part.
		for (Index k=0; k<9; ++k)
		ms[NE]=ms[SW]=0.5*(stencilL[position[SW]]+stencilL[position[NE]]);
		ms[N]=ms[S]=0.5*(stencilL[position[S]]+stencilL[position[N]]);
		ms[NW]=ms[SE]=0.5*(stencilL[position[SE]]+stencilL[position[NW]]);
		ms[E]=ms[W]=0.5*(stencilL[position[W]]+stencilL[position[E]]);
		ms[C]=stencilL[position[C]];
		symsum+=2*ms[NE]+2*ms[N]+2*ms[NW]+2*ms[E]+ms[C];
		mt[SW]=0.5*(stencilL[position[SW]]-stencilL[position[NE]]);
		mt[S]=0.5*(stencilL[position[S]]-stencilL[position[N]]);
		mt[SE]=0.5*(stencilL[position[SE]]-stencilL[position[NW]]);
		mt[W]=0.5*(stencilL[position[W]]-stencilL[position[E]]);
		mt[NE]=0.5*(stencilL[position[NE]]-stencilL[position[SW]]);
		mt[N]=0.5*(stencilL[position[N]]-stencilL[position[S]]);
		mt[NW]=0.5*(stencilL[position[NW]]-stencilL[position[SE]]);
		mt[E]=0.5*(stencilL[position[E]]-stencilL[position[E]]);
		mt[C]=0.0;
		d_w=std::max(std::fabs(ms[SW]+ms[W]+ms[NW]),
					 std::max(std::fabs(ms[SW]),
							  std::fabs(ms[NW])));
		d_e=std::max(std::fabs(ms[SE]+ms[E]+ms[NE]),
					 std::max(std::fabs(ms[SE]),
							  std::fabs(ms[NE])));
		d_n=std::max(std::fabs(ms[NW]+ms[N]+ms[NE]),
					 std::max(std::fabs(ms[NW]),
							  std::fabs(ms[NE])));
		d_s=std::max(std::fabs(ms[SW]+ms[S]+ms[SE]),
					 std::max(std::fabs(ms[SW]),
							  std::fabs(ms[SE])));
		sigma1=0.5*std::min(1.0, std::fabs(1-symsum/stencilL[0]));
		c_1=mt[SE]+mt[E]+mt[NE]-mt[SW]-mt[W]-mt[NW];
		w_w=sigma1*(1+(d_w-d_e)/(d_w+d_e)+c_1/(d_w+d_e+d_n+d_s));
		w_e=2*sigma1-w_w;
		weight1=std::min(2*sigma1, std::max(w_w, 0.0));
		t_[3]=weight1; 

		// S
		stencilL=stencil.getL(C,sx,sy-1,nx,ny);
		symsum=0;
	    
		// Divide the stencil defined by stencilL und position into a symmetric 
		// and an antisymmetric part.
		ms[NE]=ms[SW]=0.5*(stencilL[position[SW]]+stencilL[position[NE]]);
		ms[N]=ms[S]=0.5*(stencilL[position[S]]+stencilL[position[N]]);
		ms[NW]=ms[SE]=0.5*(stencilL[position[SE]]+stencilL[position[NW]]);
		ms[E]=ms[W]=0.5*(stencilL[position[W]]+stencilL[position[E]]);
		ms[C]=stencilL[position[C]];
		symsum+=2*ms[NE]+2*ms[N]+2*ms[NW]+2*ms[E]+ms[C];
		mt[SW]=0.5*(stencilL[position[SW]]-stencilL[position[NE]]);
		mt[S]=0.5*(stencilL[position[S]]-stencilL[position[N]]);
		mt[SE]=0.5*(stencilL[position[SE]]-stencilL[position[NW]]);
		mt[W]=0.5*(stencilL[position[W]]-stencilL[position[E]]);
		mt[NE]=0.5*(stencilL[position[NE]]-stencilL[position[SW]]);
		mt[N]=0.5*(stencilL[position[N]]-stencilL[position[S]]);
		mt[NW]=0.5*(stencilL[position[NW]]-stencilL[position[SE]]);
		mt[E]=0.5*(stencilL[position[E]]-stencilL[position[E]]);
		mt[C]=0.0;
		d_w=std::max(std::fabs(ms[SW]+ms[W]+ms[NW]),
					 std::max(std::fabs(ms[SW]),
							  std::fabs(ms[NW])));
		d_e=std::max(std::fabs(ms[SE]+ms[E]+ms[NE]),
					 std::max(std::fabs(ms[SE]),
							  std::fabs(ms[NE])));
		d_n=std::max(std::fabs(ms[NW]+ms[N]+ms[NE]),
					 std::max(std::fabs(ms[NW]),
							  std::fabs(ms[NE])));
		d_s=std::max(std::fabs(ms[SW]+ms[S]+ms[SE]),
					 std::max(std::fabs(ms[SW]),
							  std::fabs(ms[SE])));
		sigma2=0.5*std::min(1.0, std::fabs(1-symsum/stencilL[0]));
		c_2=mt[NW]+mt[N]+mt[NE]-mt[SW]-mt[S]-mt[SE];
		w_n=sigma2*(1+(d_s-d_n)/(d_s+d_n)+c_2/(d_w+d_e+d_n+d_s));
		w_s=2*sigma2-w_n;
		weight2=std::min(2*sigma2, std::max(w_n, 0.0));
		t_[4]=weight2;

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

	}
    const Precision t[] = {
    1.0,   1.0/2, 1.0/2,
    1.0/2, 1.0/2, 1.0/4,
    1.0/4, 1.0/4, 1.0/4};
    t_ = NumericArray(t,9);
    return t_;
}

}
