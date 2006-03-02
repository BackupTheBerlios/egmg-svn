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
    for (size_t j=0; j<=ny; ++j)
        for (size_t i=0; i<=nx; ++i)
            result[2*j*(nxNew+1)+2*i]=u[j*(nx+1)+i];

    //interpolation of fine grid points on coarse grid lines
    for (size_t j=0; j<=nyNew; j+=2)
        for (size_t i=1; i<=nxNew; i+=2)
            result[j*(nxNew+1)+i]=
                1./2*(result[j*(nxNew+1)+i-1]
                     +result[j*(nxNew+1)+i+1]);

    //interpolation of fine grid points on fine grid lines and coarse
    //grid columns
    for (size_t j=1; j<=nyNew; j+=2)
        for (size_t i=0; i<=nxNew; i+=2)
            result[j*(nxNew+1)+i]=
                1./2*(result[(j-1)*(nxNew+1)+i]
                     +result[(j+1)*(nxNew+1)+i]);

    //interpolation of fine grid points on fine grid lines and fine grid columns
    for (size_t j=1; j<=nyNew; j+=2)
        for (size_t i=1; i<=nxNew; i+=2)
            result[j*(nxNew+1)+i]=
                1./2*(result[(j-1)*(nxNew+1)+i-1]
                     +result[(j+1)*(nxNew+1)+i+1]);

    return result;
}

}
