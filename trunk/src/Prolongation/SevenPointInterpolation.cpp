/** \file SevenPointInterpolation.cpp
 * \author Benedikt Engbroks
 * \brief Contains the implementaion of the class SevenPointInterpolation
 * \see SevenPointInterpolation.h
 */
#include "SevenPointInterpolation.h"

namespace mg
{
    
NumericArray SevenPointInterpolation::prolongate(
    const NumericArray& u,
    const Stencil&,
    const size_t nx,
    const size_t ny) const
{
    const size_t nxNew = 2*nx;
    const size_t nyNew = 2*ny;
    NumericArray result((nxNew+1)*(nyNew+1));

    //"interpolation" of coarse grid points
    for (size_t sy=0; sy<=ny; ++sy)
        for (size_t sx=0; sx<=nx; ++sx)
            result[2*sy*(nxNew+1)+2*sx]=u[sy*(nx+1)+sx];

    //interpolation of fine grid points on coarse grid lines
    for (size_t sy=0; sy<=nyNew; sy+=2)
        for (size_t sx=1; sx<=nxNew; sx+=2)
            result[sy*(nxNew+1)+sx]=
                1./2*(result[sy*(nxNew+1)+sx-1]
                     +result[sy*(nxNew+1)+sx+1]);

    //interpolation of fine grid points on fine grid lines and coarse
    //grid columns
    for (size_t sy=1; sy<=nyNew; sy+=2)
        for (size_t sx=0; sx<=nxNew; sx+=2)
            result[sy*(nxNew+1)+sx]=
                1./2*(result[(sy-1)*(nxNew+1)+sx]
                     +result[(sy+1)*(nxNew+1)+sx]);

    //interpolation of fine grid points on fine grid lines and fine grid columns
    for (size_t sy=1; sy<=nyNew; sy+=2)
        for (size_t sx=1; sx<=nxNew; sx+=2)
            result[sy*(nxNew+1)+sx]=
                1./2*(result[(sy-1)*(nxNew+1)+sx-1]
                     +result[(sy+1)*(nxNew+1)+sx+1]);

    return result;
}

}
