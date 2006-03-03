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
    const size_t nx,
    const size_t ny) const
{
    register const size_t nxNew=2*nx;
    register const size_t nyNew=2*ny;
    NumericArray result((nxNew+1)*(nyNew+1));
    PositionArray jx=stencil.getJx(C);
    PositionArray jy=stencil.getJy(C);
    NumericArray stencilL=stencil.getL(C,0,0,nx,ny);
    std::valarray<size_t> position(9);
    NumericArray ms(9);
    NumericArray mt(9);
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
    for (size_t sy=0; sy<=ny; ++sy)
        for (size_t sx=0; sx<=nx; ++sx)
            result[2*sy*(nxNew+1)+2*sx]=u[sy*(nx+1)+sx];
    
    //interpolation of fine grid points on coarse grid lines
    for (size_t sy=0; sy<=nyNew; sy+=2)
        for (size_t sx=1; sx<=nxNew; sx+=2)
        {
            stencilL=stencil.getL(C,sx,sy,nx,ny);
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
            result[sy*(nxNew+1)+sx]=
                 weight1*result[sy*(nxNew+1)+sx-1]
                +weight2*result[sy*(nxNew+1)+sx+1];
        }
    
    // interpolation of fine grid points on fine grid lines and coarse grid 
    // columns
    for (size_t sy=1; sy<=nyNew; sy+=2)
        for (size_t sx=0; sx<=nxNew; sx+=2)
        {
            stencilL=stencil.getL(C,sx,sy,nx,ny);
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
            result[sy*(nxNew+1)+sx]=
                 weight1*result[(sy-1)*(nxNew+1)+sx]
                +weight2*result[(sy+1)*(nxNew+1)+sx];
        }
            
    //interpolation of fine grid points on fine grid lines and fine grid columns
    for (size_t sy=1; sy<=nyNew; sy+=2)
        for (size_t sx=1; sx<=nxNew; sx+=2)
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
const NumericArray& DeZeeuwInterpolation::get_I(
    const size_t sx,
    const size_t sy, 
    const size_t nx,
    const size_t ny,
    const Stencil& stencil)
{
    NumericArray stencilL=stencil.getL(C,0,0,nx,ny);
    PositionArray jx=stencil.getJx(C);
    PositionArray jy=stencil.getJy(C);
    std::valarray<size_t> position(9);
    for (size_t jj=0;jj<jx.size();jj++)
    {
        position[(jx[jj]+1)+3*(jy[jj]+1)]=jj;
    }
    NumericArray mS(9);
    NumericArray mT(9);

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
    stencilL=stencil.getL(C,sx-1,sy,nx,ny);
    symsum=0;
    
    // Divide the stencil defined by stencilL und position into a symmetric 
    // and an antisymmetric part.
    for (size_t k=0; k<9; ++k)
    {
        mS[k]=0.5*(stencilL[position[k]]+stencilL[position[8-k]]);
        symsum+=mS[k];
        mT[k]=0.5*(stencilL[position[k]]-stencilL[position[8-k]]);
    }
    d_w=std::max(std::fabs(mS[SW]+mS[W]+mS[NW]),
                 std::max(std::fabs(mS[SW]),
                          std::fabs(mS[NW])));
    d_e=std::max(std::fabs(mS[SE]+mS[E]+mS[NE]),
                 std::max(std::fabs(mS[SE]),
                          std::fabs(mS[NE])));
    d_n=std::max(std::fabs(mS[NW]+mS[N]+mS[NE]),
                 std::max(std::fabs(mS[NW]),
                          std::fabs(mS[NE])));
    d_s=std::max(std::fabs(mS[SW]+mS[S]+mS[SE]),
                 std::max(std::fabs(mS[SW]),
                          std::fabs(mS[SE])));
    sigma1=0.5*std::min(1.0, std::fabs(1-symsum/stencilL[0]));
    c_1=mT[SE]+mT[E]+mT[NE]-mT[SW]-mT[W]-mT[NW];
    w_w=sigma1*(1+(d_w-d_e)/(d_w+d_e)+c_1/(d_w+d_e+d_n+d_s));
    w_e=2*sigma1-w_w;
    weight2=std::min(2*sigma1, std::max(w_e, 0.0));
    t_[1]=weight2; 

    // N
    stencilL=stencil.getL(C,sx,sy+1,nx,ny);
    symsum=0;
    
    // Divide the stencil defined by stencilL und position into a symmetric 
    // and an antisymmetric part.
    for (int k=0; k<9; ++k)
    {
        mS[k]=0.5*(stencilL[position[k]]+stencilL[position[8-k]]);
        symsum+=mS[k];
        mT[k]=0.5*(stencilL[position[k]]-stencilL[position[8-k]]);
    }
    d_w=std::max(std::fabs(mS[SW]+mS[W]+mS[NW]),
                 std::max(std::fabs(mS[SW]),
                          std::fabs(mS[NW])));
    d_e=std::max(std::fabs(mS[SE]+mS[E]+mS[NE]),
                 std::max(std::fabs(mS[SE]),
                          std::fabs(mS[NE])));
    d_n=std::max(std::fabs(mS[NW]+mS[N]+mS[NE]),
                 std::max(std::fabs(mS[NW]),
                          std::fabs(mS[NE])));
    d_s=std::max(std::fabs(mS[SW]+mS[S]+mS[SE]),
                 std::max(std::fabs(mS[SW]),
                          std::fabs(mS[SE])));
    sigma2=0.5*std::min(1.0, std::fabs(1-symsum/stencilL[0]));
    c_2=mT[NW]+mT[N]+mT[NE]-mT[SW]-mT[S]-mT[SE];
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
    for (size_t k=0; k<9; ++k)
    {
        mS[k]=0.5*(stencilL[position[k]]+stencilL[position[8-k]]);
        symsum+=mS[k];
        mT[k]=0.5*(stencilL[position[k]]-stencilL[position[8-k]]);
    }
    d_w=std::max(std::fabs(mS[SW]+mS[W]+mS[NW]),
                 std::max(std::fabs(mS[SW]),
                          std::fabs(mS[NW])));
    d_e=std::max(std::fabs(mS[SE]+mS[E]+mS[NE]),
                 std::max(std::fabs(mS[SE]),
                          std::fabs(mS[NE])));
    d_n=std::max(std::fabs(mS[NW]+mS[N]+mS[NE]),
                 std::max(std::fabs(mS[NW]),
                          std::fabs(mS[NE])));
    d_s=std::max(std::fabs(mS[SW]+mS[S]+mS[SE]),
                 std::max(std::fabs(mS[SW]),
                          std::fabs(mS[SE])));
    sigma1=0.5*std::min(1.0, std::fabs(1-symsum/stencilL[0]));
    c_1=mT[SE]+mT[E]+mT[NE]-mT[SW]-mT[W]-mT[NW];
    w_w=sigma1*(1+(d_w-d_e)/(d_w+d_e)+c_1/(d_w+d_e+d_n+d_s));
    w_e=2*sigma1-w_w;
    weight1=std::min(2*sigma1, std::max(w_w, 0.0));
    t_[3]=weight1; 

    // S
    stencilL=stencil.getL(C,sx,sy-1,nx,ny);
    symsum=0;
    
    // Divide the stencil defined by stencilL und position into a symmetric 
    // and an antisymmetric part.
    for (int k=0; k<9; ++k)
    {
        mS[k]=0.5*(stencilL[position[k]]+stencilL[position[8-k]]);
        symsum+=mS[k];
        mT[k]=0.5*(stencilL[position[k]]-stencilL[position[8-k]]);
    }
    d_w=std::max(std::fabs(mS[SW]+mS[W]+mS[NW]),
                 std::max(std::fabs(mS[SW]),
                          std::fabs(mS[NW])));
    d_e=std::max(std::fabs(mS[SE]+mS[E]+mS[NE]),
                 std::max(std::fabs(mS[SE]),
                          std::fabs(mS[NE])));
    d_n=std::max(std::fabs(mS[NW]+mS[N]+mS[NE]),
                 std::max(std::fabs(mS[NW]),
                          std::fabs(mS[NE])));
    d_s=std::max(std::fabs(mS[SW]+mS[S]+mS[SE]),
                 std::max(std::fabs(mS[SW]),
                          std::fabs(mS[SE])));
    sigma2=0.5*std::min(1.0, std::fabs(1-symsum/stencilL[0]));
    c_2=mT[NW]+mT[N]+mT[NE]-mT[SW]-mT[S]-mT[SE];
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

    return t_;   
}

}
