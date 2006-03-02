/** \file DendyInterpolation.cpp
 * \author Benedikt Engbroks
 * \brief Contains the implementation of the class DendyInterpolation
 * \see DendyInterpolation.h
 */
#include "DendyInterpolation.h"

namespace mg
{
// This function only works for max. compact 9-point stencils
NumericArray DendyInterpolation::prolongate(
    const NumericArray& u,
    const Stencil& stencil,
    const size_t nx,
    const size_t ny) const
{
    const size_t nxNew=2*nx;
    const size_t nyNew=2*ny;
    NumericArray result((nxNew+1)*(nyNew+1));
    PositionArray jx=stencil.getJx(C);
    PositionArray jy=stencil.getJy(C);
    NumericArray stencilL=NumericArray();
    std::valarray<size_t> position(9);
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
        position[(jx[j]+1)+3*(jy[j]+1)]=j;

    Precision scale=0;
    Precision weight1=0;
    Precision weight2=0;
    Precision erg=0;
    
    //"interpolation" of coarse grid points
    for (size_t j=0; j<=ny; ++j)
        for (size_t i=0; i<=nx; ++i)
            result[2*j*(nxNew+1)+2*i]=u[j*(nx+1)+i];
    
    //interpolation of fine grid points on coarse grid lines
    for (size_t j=0; j<=nyNew; j+=2)
        for (size_t i=1; i<=nxNew; i+=2)
        {
            stencilL=stencil.getL(C,i,j,nx,ny);
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
            result[j*(nxNew+1)+i]=(weight1*result[j*(nxNew+1)+i-1]+weight2*result[j*(nxNew+1)+i+1])/scale;
        }
    
    //interpolation of fine grid points on fine grid lines and coarse grid columns
    for (size_t j=1; j<=nyNew; j+=2)
        for (size_t i=0; i<=nxNew; i+=2)
        {
            stencilL=stencil.getL(C,i,j,nx,ny);
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
            result[j*(nxNew+1)+i]=(weight1*result[(j-1)*(nxNew+1)+i]+weight2*result[(j+1)*(nxNew+1)+i])/scale;
        }
            
    //interpolation of fine grid points on fine grid lines and fine grid columns
    for (size_t j=1; j<=nyNew; j+=2)
        for (size_t i=1; i<=nxNew; i+=2)
        {
            stencilL=stencil.getL(C,i,j,nx,ny);
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
