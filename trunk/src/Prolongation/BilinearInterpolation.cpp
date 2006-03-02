/** \file BilinearInterpolation.cpp
 * \author <a href="mailto:mail@jirikraus.de">Jiri Kraus</a>
 * \brief Contains the implementation of the class BilinearInterpolation
 */
#include "BilinearInterpolation.h"

namespace mg
{
    
NumericArray BilinearInterpolation::prolongate(
    const NumericArray& u,
    const Stencil&,
    const size_t nx,
    const size_t ny) const
{
    const size_t nxNew = 2*nx;
    const size_t nyNew = 2*ny;
    NumericArray result(0.0,(nxNew+1)*(nxNew+1));

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
    //interpolation of fine grid points on fine grid lines
    for (size_t j=1; j<=nyNew; j+=2)
        for (size_t i=0; i<=nxNew; ++i)
            result[j*(nxNew+1)+i]=
                1./2*(result[(j-1)*(nxNew+1)+i]
                     +result[(j+1)*(nxNew+1)+i]);
    return result;
}

}
