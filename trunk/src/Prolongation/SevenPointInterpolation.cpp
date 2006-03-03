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
    const Index nx,
    const Index ny) const
{
    const Index nxNew = 2*nx;
    const Index nyNew = 2*ny;
    NumericArray result((nxNew+1)*(nyNew+1));

    //"interpolation" of coarse grid points
    for (Index sy=0; sy<=ny; ++sy)
        for (Index sx=0; sx<=nx; ++sx)
            result[2*sy*(nxNew+1)+2*sx]=u[sy*(nx+1)+sx];

    //interpolation of fine grid points on coarse grid lines
    for (Index sy=0; sy<=nyNew; sy+=2)
        for (Index sx=1; sx<=nxNew; sx+=2)
            result[sy*(nxNew+1)+sx]=
                1./2*(result[sy*(nxNew+1)+sx-1]
                     +result[sy*(nxNew+1)+sx+1]);

    //interpolation of fine grid points on fine grid lines and coarse
    //grid columns
    for (Index sy=1; sy<=nyNew; sy+=2)
        for (Index sx=0; sx<=nxNew; sx+=2)
            result[sy*(nxNew+1)+sx]=
                1./2*(result[(sy-1)*(nxNew+1)+sx]
                     +result[(sy+1)*(nxNew+1)+sx]);

    //interpolation of fine grid points on fine grid lines and fine grid columns
    for (Index sy=1; sy<=nyNew; sy+=2)
        for (Index sx=1; sx<=nxNew; sx+=2)
            result[sy*(nxNew+1)+sx]=
                1./2*(result[(sy-1)*(nxNew+1)+sx-1]
                     +result[(sy+1)*(nxNew+1)+sx+1]);

    return result;
}

}
