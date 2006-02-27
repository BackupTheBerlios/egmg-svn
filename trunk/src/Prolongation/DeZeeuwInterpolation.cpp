/** \file DeZeeuwInterpolation.cpp
 * \author Benedikt Engbroks
 * \brief Contains the implementation of the class DeZeeuwInterpolation
 * \see DeZeeuwInterpolation.h
 */

//TODO: This file needs cleanup

#include "DeZeeuwInterpolation.h"
#include "../Stencil/Stencil.h"
#include <algorithm>

namespace mg
{
// This function only works for max. compact 9-point stencils
std::valarray<Precision> DeZeeuwInterpolation::prolongate(
    const std::valarray<Precision>& u,
    const Stencil& stencil,
    const size_t nx,
    const size_t ny) const
{
    register const size_t nxNew=2*nx;
    register const size_t nyNew=2*ny;
    std::valarray<Precision> result((nxNew+1)*(nyNew+1));
    std::valarray<int> jx=stencil.getJx(c);
    std::valarray<int> jy=stencil.getJy(c);
    std::valarray<Precision> stencilL=stencil.getL(c,0,0,nx,ny);
    std::valarray<size_t> position(9);
    std::valarray<Precision> ms(9);
    std::valarray<Precision> mt(9);
    // position[0]=No. of the sw-element of the stencil
    // position[1]=No. of the s-element of the stencil
    // position[2]=No. of the se-element of the stencil
    // position[3]=No. of the w-element of the stencil
    // position[4]=No. of the c-element of the stencil (i.e. 0)
    // position[5]=No. of the e-element of the stencil
    // position[6]=No. of the nw-element of the stencil
    // position[7]=No. of the n-element of the stencil
    // position[8]=No. of the ne-element of the stencil
    // (see also the definitions of the constants in the .h-file)
    for (size_t j=0; j<jx.size(); ++j)
    {
        position[(jx[j]+1)+3*(jy[j]+1)]=j;
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

    //"interpolation" of coarse grid points
    for (size_t j=0; j<=ny; ++j)
        for (size_t i=0; i<=nx; ++i)
            result[2*j*(nxNew+1)+2*i]=u[j*(nx+1)+i];
    
    //interpolation of fine grid points on coarse grid lines
    for (size_t j=0; j<=nyNew; j+=2)
        for (size_t i=1; i<=nxNew; i+=2)
        {
            stencilL=stencil.getL(c,i,j,nx,ny);
            symsum=0;
            
            // Divide the stencil defined by stencilL und position into a
            // symmetric and an antisymmetric part.
            for (size_t k=0; k<9; ++k)
            {
                ms[k]=0.5*(stencilL[position[k]]+stencilL[position[8-k]]);
                symsum+=ms[k];
                mt[k]=0.5*(stencilL[position[k]]-stencilL[position[8-k]]);
            }
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
            result[j*(nxNew+1)+i]=
                 weight1*result[j*(nxNew+1)+i-1]
                +weight2*result[j*(nxNew+1)+i+1];
        }
    
    // interpolation of fine grid points on fine grid lines and coarse grid 
    // columns
    for (size_t j=1; j<=nyNew; j+=2)
        for (size_t i=0; i<=nxNew; i+=2)
        {
            stencilL=stencil.getL(c,i,j,nx,ny);
            symsum=0;
            
            // Divide the stencil defined by stencilL und position into a
            // symmetric and an antisymmetric part.
            for (int k=0; k<9; ++k)
            {
                ms[k]=0.5*(stencilL[position[k]]+stencilL[position[8-k]]);
                symsum+=ms[k];
                mt[k]=0.5*(stencilL[position[k]]-stencilL[position[8-k]]);
            }
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
            result[j*(nxNew+1)+i]=
                 weight1*result[(j-1)*(nxNew+1)+i]
                +weight2*result[(j+1)*(nxNew+1)+i];
        }
            
    //interpolation of fine grid points on fine grid lines and fine grid columns
    for (size_t j=1; j<=nyNew; j+=2)
        for (size_t i=1; i<=nxNew; i+=2)
        {
            stencilL=stencil.getL(c,i,j,nx,ny);
            erg=0;
            scale=-stencilL[0];
            if (position[W]!=0)
                erg+=stencilL[position[W]]*result[j*(nxNew+1)+i-1];
            if (position[E]!=0)
                erg+=stencilL[position[E]]*result[j*(nxNew+1)+i+1];
            if (position[S]!=0)
                erg+=stencilL[position[S]]*result[(j-1)*(nxNew+1)+i];
            if (position[SW]!=0)
                erg+=stencilL[position[SW]]*result[(j-1)*(nxNew+1)+i-1];
            if (position[SE]!=0)
                erg+=stencilL[position[SE]]*result[(j-1)*(nxNew+1)+i+1];
            if (position[N]!=0)
                erg+=stencilL[position[N]]*result[(j+1)*(nxNew+1)+i];
            if (position[NW]!=0)
                erg+=stencilL[position[NW]]*result[(j+1)*(nxNew+1)+i-1];
            if (position[NE]!=0)
                erg+=stencilL[position[NE]]*result[(j+1)*(nxNew+1)+i+1];
            result[j*(nxNew+1)+i]=erg/scale;
        }
    return result;
}


}
